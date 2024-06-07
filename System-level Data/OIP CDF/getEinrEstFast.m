% Obtain EINR estimation for each range-velocity bins
function fast_einr_est = getEinrEstFast(Y_rD_3D,nCenter,lCenter,param,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells)
% parameters 
n_min = max(nCenter-rangeGI-rangeTrainCells,1);
n_max = min(nCenter+rangeGI+rangeTrainCells,param.Nrfft);
l_min = max(lCenter-velocityGI-velocityTrainCells,1);
l_max = min(lCenter+velocityGI+velocityTrainCells,param.Nvfft);
numNearbyRangebins = n_max-n_min+1;
numNearbyVelocitybins = l_max-l_min+1;
noiseEstRaw = zeros(numNearbyRangebins,numNearbyVelocitybins); % esimation of noise power
b_tuta_estRaw = zeros(numNearbyRangebins,numNearbyVelocitybins,param.Nt,param.numInt); % esimation of b tuta on different range-velocity bins 
hsquare_estRaw = zeros(numNearbyRangebins,numNearbyVelocitybins,param.Nt,param.numInt);% esimation of h sqaure on different range-velocity bins
% GS EINR estimation
for n = n_min:n_max
    for l = l_min:l_max
        Y_rD_3D_n_l = Y_rD_3D(n,:,l)';
        % Obtain projected signal y_1, y_2, ..., y_M 
        y_rD_3D_n_l = reshape(Y_rD_3D_n_l,[param.Nr,param.Nt]); % decomponse y into y_1, y_2, ..., y_M
        Y_rD_3D_n_l_proj = reshape(param.P_Ar_int_orth*y_rD_3D_n_l,[],1);
        noiseEstRaw(n-n_min+1,l-l_min+1) = 2*norm(Y_rD_3D_n_l_proj)^2/(param.Nt*(param.Nr-param.numInt)); % estimator of noise = 2*||[P_Ar_orth*y_1;.. P_Ar_orth*y_M]||^2/(MN-MQ)
        for q = 1:param.numInt
            b_q = param.B_RS(:,q); % the q-th column of B
            at_int_est_n_l_q = reshape((param.Ar_int*b_q)'*y_rD_3D_n_l,[],1); % estimator of at_int_q = [(param.Ar_int*b_q)'*y_1; ...; (param.Ar_int*b_q)'*y_M]
            b_tuta_estRaw(n-n_min+1,l-l_min+1,:,q) = fft(at_int_est_n_l_q,param.Nt)/param.Nt;
        end
        hsquare_estRaw(n-n_min+1,l-l_min+1,:,:) = abs(squeeze(b_tuta_estRaw(n-n_min+1,l-l_min+1,:,:))).^2;
    end
end
% EINR estimation
fast_einr_est= getEinrEstSingleBin(nCenter-n_min+1,lCenter-l_min+1,hsquare_estRaw,noiseEstRaw,param,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells);
end
