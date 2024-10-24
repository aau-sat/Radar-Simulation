function [time_data_sim, time_data_det] = generate_3D_timedomain_DPN(intensity,distance, velocity, c0, samplesNumber, chirpsNumber,nvah,nvav,fs, LowerChirpFrequency, ChirpRate, az_vec, el_vec, spacing, delta_t,pn_mask_gain, pn_mask_freq)
    % Init variables:
    lambda = c0/LowerChirpFrequency; 
    ptot = 1:chirpsNumber;
    %PRT = 1.0192e-04; 
    % Tc = PRT;
    Tc = samplesNumber/fs; 
    n_tgs = length(distance);
    
    df = 1/Tc;
    f_vec_PSD = df:df:fs/2;
    noise_len = length(f_vec_PSD);
    
    
    tau = 2*(distance + velocity.*delta_t)./c0 + ...
        LowerChirpFrequency*spacing*sin(az_vec)/c0 + ...
        LowerChirpFrequency*spacing*sin(el_vec)/c0;
    
    % Generate DPN (not optimized)
    for a = 1:(nvah + nvav)
       Noise_vec = 1*ones(1,noise_len); % adjust
       Noise_vec = Noise_vec.* exp(1i*2*pi*rand(1,noise_len)); 
       for tt = 1:n_tgs
          [dpn_gen] = generate_dpn_samples_mod(pn_mask_gain, pn_mask_freq, f_vec_PSD, tau(tt), Noise_vec); % 1 x Ns
           for p = 1:chirpsNumber
           dpn(tt,p,:,a) = dpn_gen;
           end
       end 
    end
    
    
    
    % Signal gen
    t = [1:samplesNumber]/fs;
    
    fd = 2*velocity./lambda;
    
    % all matrices:
    f1a = (2*ChirpRate*distance./c0 + fd);  % here need the kronecker
    f1b = kron(t,f1a);
    f1 = reshape(f1b,[length(distance),1,samplesNumber,1]);
    
    f2a = LowerChirpFrequency*spacing*sin(az_vec)/c0;
    f2b = kron([1:1:nvah]-1,f2a); % note that we consider antenna 0 as the first (not 1)
    f2c = LowerChirpFrequency*spacing*sin(el_vec)/c0;
    f2d = kron([1:1:nvav]-1,f2c); % note that we consider antenna 0 as the first (not 1)
    f2tot = cat(2,f2b,f2d);
    f2 = reshape(f2tot,[n_tgs,1,1,nvah+nvav]); % contribution over all antennas
    
    f3a = 2*fd*Tc;
    f3b = kron(ptot,f3a);
    f3 = reshape(f3b,[length(distance),chirpsNumber,1,1]); 

    f4 = (2*LowerChirpFrequency*distance/c0)';

    costerm = 2*pi*(f1+f2+f3+f4);
    beta = 1; % adjust according to needed intensity/amplitude
    x_if_s = beta*intensity.*cos(costerm + dpn);
    x_if_s_det = beta*intensity.*cos(costerm); 
    time_data_sim = (sum(x_if_s,1,'omitnan'));   
    time_data_det = (sum(x_if_s_det,1,'omitnan')); 
end


