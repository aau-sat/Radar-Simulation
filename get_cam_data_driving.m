function [seg_clouds,xpc_mat,ypc_mat,zpc_mat,amplitude_ToF,amplitude_radar,cloud_toSave] = get_cam_data_driving(tcpipClient,NrPixel_x, NrPixel_y, fov_horizontal,fov_vertical )
 
    % Read data from Unity
    fwrite(tcpipClient,'READY'); % send ready command to Unity
    image_length = fread(tcpipClient, 4); %read length of image to be send
    len = (typecast(uint8(image_length)','int32')); 
    
    % Comment or uncomment according to version
    % utf8Bytes = fread(tcpipClient, 46); %read length of the timestpamp
    % timeStampStr=convertCharsToStrings(native2unicode(utf8Bytes, 'UTF-8')); %#ok<N2UNI>

    rawData = fread(tcpipClient, double(len)); % read image bytes in png format
    rawData_uint = uint8(rawData);                      % cast them to "uint8" if they are not already
    rawData_float = typecast( fliplr(rawData_uint) , 'single');
    raw_Data_double = double(rawData_float);

    % color channels:
    R = raw_Data_double(1:3:end-1);
    R = reshape(R, NrPixel_x, NrPixel_y);
    R = R';
    amplitude_ToF = R(end:-1:1,:); % RED - previously used for ToF amplitude  
    G = raw_Data_double(2:3:end-1);
    G = reshape(G, NrPixel_x, NrPixel_y);
    G = G';
    depth = G(end:-1:1,:); % GREEN - depth     
    B = raw_Data_double(3:3:end);
    B = reshape(B, NrPixel_x, NrPixel_y);
    B = B';
    amplitude_radar = B(end:-1:1,:);% BLUE - Radar amplitude

    %time_stamp = raw_Data_double(end); % Not used

    % Camera coordinates computation (i.e., GT point cloud -> sort of perfect LiDaR)
    X_max_SS = tan(deg2rad(fov_horizontal/2)); % max X
    Y_max_SS = tan(deg2rad(fov_vertical/2)); % max Y
    Z = depth'; 
    X = (Z*X_max_SS) .* (((repmat([1:NrPixel_x]',1,NrPixel_y)) - (NrPixel_x/2))/(NrPixel_x/2));
    Y = (Z*Y_max_SS) .* (((repmat([1:NrPixel_y],NrPixel_x,1)) - (NrPixel_y/2))/(NrPixel_y/2));
    GT_pc(:,:,1) = X;
    GT_pc(:,:,2) = Y;
    GT_pc(:,:,3) = Z;
    GT_pc(:,:,4) = 1; % if needed for transforms

    % Coordinate change (Note, this now will have Z up and Y forward):
    xpc_mat = GT_pc(:,:,1);
    ypc_mat = GT_pc(:,:,3);
    zpc_mat = -GT_pc(:,:,2);
 
    % Pre-select points based on amplitude. Too low intensity are from
    % targets that anyway will never be detected
    int_thr = 0.3; % adjust
    [idxs] = find(amplitude_radar'>=int_thr); 
    xpc_mat = xpc_mat(idxs);
    ypc_mat = ypc_mat(idxs);
    zpc_mat = zpc_mat(idxs);
    int_mat = amplitude_radar(idxs);

    % Also remove unwanted points based on specific experiment, if needed
    % zpc_mat(zpc_mat<-1.2) = NaN;
    % zpc_mat(zpc_mat>2.25) = NaN;
    % xpc_mat(abs(xpc_mat)>5.5) = NaN;
    
    cloud_pos = pointCloud(cat(3,xpc_mat, ypc_mat,zpc_mat)); % create RGB point cloud
    cloud_pos.Intensity = int_mat;
    cloud_pos = removeInvalidPoints(cloud_pos); % removes nans
    
    % first downsample
    ds_pc = pcdownsample(cloud_pos,'gridAverage',0.5,'PreserveStructure',true); %0.5 0.25
    cloud_toSave = ds_pc;

   
    %Cluster the point cloud with a minimum of 20 points per cluster.
    minDistance = 0.75; 
    minPoints = 4; 
    [labels,numClusters] = pcsegdist(ds_pc,minDistance,'NumClusterPoints',minPoints);

    %Remove the points with a label value of 0.
    % idxValidPoints = find(labels);
    % labelColorIndex = labels(idxValidPoints);
    % segmentedPtCloud = select(ds_pc,idxValidPoints);
    for n = 1:numClusters
        seg_clouds{n} = select(ds_pc,find(labels==n));
    end 

    % figure(99)
    % clf
    % colormap(hsv(numClusters))
    % pcshow(segmentedPtCloud.Location,labelColorIndex, 'MarkerSize', 70)
    % title('Point Cloud Clusters')

end

