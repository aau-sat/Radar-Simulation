% Code used for the 3D driving scenario in the parking garage simulation
% (tested with the 032023 project and 2D CFAR)
% The code connects to Unity via TCPIP and provides:
% - time domain data with noise
% - range FFT
% - doppler FFT
% - Range Doppler Map
% - Range Angle Map
% - 3D radar point cloud

clear all, close all

% UNITY VARIABLES
NrPixel_x = 300; 
NrPixel_y = 171;
fov_vertical = 68; 
fov_horizontal = 100; 
bytes_to_receive =  3 * 4 * NrPixel_x * NrPixel_y; % max bytes to receive in case of two 32-bit color channels
tcpipClient = tcpclient('127.0.0.1',8053);
fclose(tcpipClient);
set(tcpipClient,'Timeout',10);
tcpipClient.InputBufferSize = bytes_to_receive;
fopen(tcpipClient);

% ROS SETTINGS
rosshutdown
rosinit
data_vel = rossubscriber("/vel"); % to get velocity data
pause(1);

% RADAR PARAMETERS
c0 = 300e6;
samplesNumber = 512;
chirpsNumber = 64; 
fs = 6.25e6; % ADC sampling frequency - adjust according to device 6.25 - 3.36
bandwidth =  9.03e8; % - adjust according to device 1.5e9 - 3.7
LowerChirpFrequency = 77.4e9; % - adjust according to device 76 - 77
ChirpRate = (bandwidth *fs) /samplesNumber;
Tc = samplesNumber/fs; % Chirp time
PRT = 1.0192e-04; % Pulse Repetition Rate - for max velocity calculation
max_vel = lambda/(4*PRT);
K = bandwidth/Tc; % slope of a chirp (chirp rate)
lambda = c0/LowerChirpFrequency; % wavelength
spacing = 0.6*lambda;  % spacing between receivers - adjust according to device
% Range-Doppler-Angle Axes:
xscale  = 2*K /c0; 
Zeropad = 2^(nextpow2(samplesNumber)); 
range_axis  = linspace(0,1-1/Zeropad,Zeropad)*fs/xscale;
Zeropad_velocity = Zeropad/2;
vel_axis = linspace(-max_vel/2, max_vel/2, Zeropad_velocity + 1); 
Nangle_fft = length(vel_axis);
max_angle = asind(lambda/((spacing)*2));  
angles_axis = linspace(-max_angle, max_angle, Nangle_fft);

ulah = phased.ULA('NumElements',12,'ElementSpacing',spacing);
music_est_h = phased.MUSICEstimator('SensorArray',ulah,...
            'OperatingFrequency',LowerChirpFrequency,'ScanAngles',angles_axis,...
            'DOAOutputPort',true,'NumSignalsSource','Property',...
        'NumSignals',1);
ulav = phased.ULA('NumElements',8,'ElementSpacing',spacing);
music_est_v  = phased.MUSICEstimator('SensorArray',ulav,...
        'OperatingFrequency',LowerChirpFrequency,'ScanAngles',angles_axis,...
        'DOAOutputPort',true,'NumSignalsSource','Property',...
    'NumSignals',1);


% Other Parameters used for processing:
[b,a] = butter(4,0.06,'high'); % DC comp filter, if needed
Hann_window_tr = hann(samplesNumber,'periodic');                                                  
ScaleHannWin_tr = 1/sum(Hann_window_tr);  
Cheb_window_tr = chebwin(chirpsNumber ,80);                                               
ScaleChebWin_tr = 1/sum(Cheb_window_tr ); 

% pn_mask_freq = []; % Phase Noise data are private: please contact maintainers if needed
% pn_mask_gain =  [];% Phase Noise data are private: please contact maintainers if needed

%% Processing loop
tStart = tic;
for t = 1:1000 % adjust
    tTot = tic;

    [seg_clouds,~,~,~,~,amplitude_radar, cloud_toSave] = get_cam_data_driving(tcpipClient,NrPixel_x,NrPixel_y,fov_vertical,fov_horizontal);
    seg_clouds_array = []; % array of segmented, downsampled point clouds (will need to be concatenated)

    % Select the strongest N points out of each cluster, with N = 25% of
    % all points (adjust if needed)
    for n = 1:length(seg_clouds)
        [sorted_int, idx_sort] = sort(seg_clouds{n}.Intensity, 'descend');
        seg_clouds_highest_int = select(seg_clouds{n},idx_sort(1:round(0.25*length(idx_sort))));
        seg_clouds_arr = [seg_clouds_arr, seg_clouds_highest_int];
    end

    % get linear velocity data
    velocity_msg = receive(data_vel,1);
    velocity_camera = [velocity_msg.Linear.X; velocity_msg.Linear.Y; velocity_msg.Linear.Z];   
    % then I know every 3d point for each point in the lidar cloud has
    % the negative this velocity
    v_pc = velocity_camera;

    % Build downsampled point cloud
    P_DS = pccat(seg_clouds_arr);

    % Rotate frame (if needed).
    % If no transformation is needed, this plots with Z up. Careful when
    % labeling the axes
    x_ros = P_DS.Location(:,1);
    y_ros = P_DS.Location(:,2);
    z_ros = P_DS.Location(:,3);
    pcros = pointCloud(cat(3,[x_ros,y_ros,z_ros])); % with transformed coordinate

    % Transformation needed for cartesian to spherical for radar frame
    xmatlab = P_DS.Location(:,2);
    ymatlab = -P_DS.Location(:,1);
    zmatlab = P_DS.Location(:,3);
    az_mat = atan2(ymatlab,xmatlab);
    el_mat = atan2(zmatlab,sqrt(xmatlab.^2 + ymatlab.^2)); 

    % Build vectors needed for time domain data gen
    dist_vec = sqrt(x_ros.^2 + y_ros.^2 + z_ros.^2);
    dist_vec(dist_vec < 0.1) = 0; % nearclip from Unity
    idxs = find(dist_vec~=0); % we are only interested in these pixels
    intens_radar = P_DS.Intensity;
  
    % Needed for velocity computation with scalar product:
    positions = [x_ros, y_ros, z_ros];
    positions = positions(idxs,:);
    normalizedPositions = positions ./ vecnorm(positions,2,2);

    intens_vec = intens_radar(idxs);
    dist_vec = dist_vec(idxs)'; 
    vel_vec = dot(repmat(v_pc,[1,size(positions,1)]), normalizedPositions');
    az_vec = az_mat(idxs);
    el_vec = el_mat(idxs);

    % Time parameters, if needed:
    delta_t = toc(tStart);
    tStart = tic;

    n_antennas_hori = 12;
    n_antennas_vert = 4;

    % Generate time domain data (with/without phase noise)
    [time_data_sim_dpn,~] = generate_3D_timedomain_DPN(intens_vec, dist_vec, vel_vec, c0, samplesNumber, ...
        chirpsNumber, n_antennas_hori,n_antennas_vert,fs, LowerChirpFrequency, ChirpRate, ...
        az_vec, el_vec, spacing,delta_t, pn_mask_gain, pn_mask_freq); % PN data here

    s_IF_mat_sum = squeeze(time_data_sim_dpn);

    % Add Thermal Noise and compute Range FFT spectrum
    snr_dB = 9;
    for a = 1:n_antennas_hori+n_antennas_vert
        % If DC compensation is needed, add here
        signal_power = rms(s_IF_mat_sum(1,:,a)).^2;  
        noise_power = signal_power / (10^(snr_dB / 10));  

        s_IF_mat_sum(:,:,a) = s_IF_mat_sum(:,:,a) + sqrt(noise_power)*randn(chirpsNumber,size(s_IF_mat_sum,2));

        % Range-FFT:
        S_IF_mat(:,:,a) = fft(s_IF_mat_sum(:,:,a).*Hann_window_tr',Zeropad,2)/samplesNumber;  
     end  

    % Process Range-FFT spectrum
    S_IF_mat_SSB_2 = abs(S_IF_mat(:,1:Zeropad/2,:))*1/sum(Hann_window_tr);  
    S_IF_mat_mean_SSB_2 = mean(S_IF_mat_SSB_2,1);  % mean over chirps
    S_IF_mat_mean_SSB_2 = mean(S_IF_mat_mean_SSB_2,3); % mean antennas

    S_IF_mat_mean_SSB_2 = movmean(S_IF_mat_mean_SSB_2,2); % LPF, if needed


    % Doppler Processing (2D FFT)
    for a = 1:n_antennas_hori+n_antennas_vert
        radar_cube(:,:,a) = fftshift(fft(S_IF_mat(:,1:end,a).*Cheb_window_tr, Zeropad_velocity+1,1),1);
        doppler_FFT_data{t}(:,:,a) = radar_cube(:,1:end/2,a);
    end
    doppler_FFT_mags = abs(doppler_FFT_data{t}(:,:,:).^2)*2/sum(Cheb_window_tr);  
    doppler_FFT_meanantennas = mean(doppler_FFT_mags,3);

    % figure(1)
    % clf
    % surf(range_axis(1:end/2),vel_axis, doppler_FFT_meanantennas, 'EdgeColor','none');
    % colormap jet
    % view(2)

    % Angle FFT for Range Angle Map
    angle_fft = fftshift(fft(S_IF_mat, Zeropad_velocity+1, 3), 3);
    angle_fft = squeeze(mean(angle_fft(:,1:end/2,:),1));
    RAM = abs(angle_fft);
    % figure(2)
    % clf
    % surf(angles_axis,range_axis(1:end/2), RAM, 'EdgeColor','none');
    % view(2)


    % 2D CFAR and AOA
    nguard = 8;
    ntrain = 12;
    cfar2D = phased.CFARDetector2D('GuardBandSize',nguard,'TrainingBandSize',ntrain,...
      'ProbabilityFalseAlarm',1e-5);
    rangeIndx = [1+(nguard+ntrain+1):length(range_axis)/2-(nguard+ntrain+1)];
    dopplerIndx = [1+(nguard+ntrain+1):length(vel_axis)-(nguard+ntrain+1)];
    [columnInds,rowInds] = meshgrid(rangeIndx(1):rangeIndx(end), dopplerIndx(1):dopplerIndx(end));
    CUTIdx = [rowInds(:) columnInds(:)]';

    % Perform 2D CFAR detection on RDM (mean of antennas):
    det_tmp = cfar2D(doppler_FFT_meanantennas,CUTIdx);

    % Detection label matrix:
    detections = zeros(size(doppler_FFT_meanantennas),'like',doppler_FFT_meanantennas);
    detections(dopplerIndx(1):dopplerIndx(end),rangeIndx(1):rangeIndx(end)) = reshape(det_tmp,[dopplerIndx(end) - dopplerIndx(1)+1,rangeIndx(end) - rangeIndx(1)+1]);

    % Extract indeces
    [rows_targets_idx,col_targets_idx] = find(detections==1);
    n_tgs = sum(detections(:)==1); % Number of points extracted from RDM


    % Set up AoA estimation with MUSIC (can also be done with FFT or
    % other methods)
    x_R = []; y_R = []; z_R = []; % cartesian coordinates of radar point cloud P_R
    angles_for_target_fft = [];
    angles_for_target_fft_v = [];
    virtual_arr_h = squeeze(cat(3,doppler_FFT_data{t}(:,:,1:12))); % Horizontal virtual array
    virtual_arr_v = squeeze(cat(3,doppler_FFT_data{t}(:,:,13:end))); % Vertical virtual array

    for kk = 1:n_tgs
        % Elevation

        % dp1 = squeeze(virtual_arr_v(rows_targets_idx(kk),col_targets_idx(kk),:));
        % ph = [unwrap(angle(dp1))];
        % ph_diff = mean(diff(ph));
        % el = asin(lambda*ph_diff/(2*pi*spacing));
        % z_R(end+1) = range_axis(col_targets_idx(kk))*sin(el);
    
        % Choose how many angles to allow (1 or more) for each target
        angles_musicv = music_est_v(squeeze(virtual_arr_v(rows_targets_idx(kk),col_targets_idx(kk),:))');
        % [~,aoa_locs_music] = findpeaks(squeeze(angles_music),'MinPeakHeight',0.75*max(squeeze(angles_music)));
        [~,aoa_locs_music_v] = max(squeeze(angles_musicv));
        for j = 1:length(aoa_locs_music_v) 
            angles_for_target_fft_v(j) = angles_axis(aoa_locs_music_v(j));
            el = -deg2rad(angles_for_target_fft_v(j)); % elevation (check sign!)
            z_R(end+1) = range_axis(col_targets_idx(kk))*sin(el); 
        end
        
        % Azimuth
        angles_music = music_est_h(squeeze(virtual_arr_h(rows_targets_idx(kk),col_targets_idx(kk),:))');
        % [aoa_pks_music,aoa_locs_music] = findpeaks(squeeze(angles_music),'MinPeakHeight',0.95*max(squeeze(angles_music)));
        [~,aoa_locs_music] = max(squeeze(angles_music));
        for j = 1:length(aoa_locs_music) 
            angles_for_target_fft(j) = angles_axis(aoa_locs_music(j));
            x_R(end+1) = range_axis(col_targets_idx(kk))*sin(deg2rad(angles_for_target_fft(j)))*cos(el);     % need to change sign here       
            y_R(end+1) = range_axis(col_targets_idx(kk))*cos(deg2rad(angles_for_target_fft(j)))*cos(el);
        end
    end

    % figure(3)
    % clf
    % P_R = pointCloud(cat(3,x_R,y_R,z_R));
    % h = pcshow(P_R,'MarkerSize',200,'ViewPlane', 'XY')
    % % hold on, pcshow(pcros,'MarkerSize',200,'ViewPlane', 'XY')
    % % xlim([-7,7]), ylim([0,40])

    toc(tTot)
end
