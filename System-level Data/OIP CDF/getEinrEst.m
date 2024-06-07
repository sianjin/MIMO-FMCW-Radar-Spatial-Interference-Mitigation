% Obtain EINR estimation for each range-velocity bins
function einr_est = getEinrEst(Y_rD_3D,param)
noiseEstRaw = zeros(param.Nrfft,param.Nvfft); % esimation of noise power
b_tuta_estRaw = zeros(param.Nrfft,param.Nvfft,param.Nt,param.numInt); % esimation of b tuta on different range-velocity bins 
hsquare_estRaw = zeros(param.Nrfft,param.Nvfft,param.Nt,param.numInt);% esimation of h sqaure on different range-velocity bins
einr_estRaw = zeros(param.Nrfft,param.Nvfft,param.Nt,param.numInt); % esimation of EINR on different range-velocity bins 
einr_est = zeros(param.Nrfft,param.Nvfft,param.Nt,param.numInt); % esimation of EINR from nearby range-velocity bins 
% GS EINR estimation
for n = 1: param.Nrfft
    for l = 1: param.Nvfft
        Y_rD_3D_n_l = Y_rD_3D(n,:,l)';
        % Obtain projected signal y_1, y_2, ..., y_M 
        y_rD_3D_n_l = reshape(Y_rD_3D_n_l,[param.Nr,param.Nt]); % decomponse y into y_1, y_2, ..., y_M
        Y_rD_3D_n_l_proj = reshape(param.P_Ar_int_orth*y_rD_3D_n_l,[],1);
        noiseEstRaw(n,l) = 2*norm(Y_rD_3D_n_l_proj)^2/(param.Nt*(param.Nr-param.numInt)); % estimator of noise = 2*||[P_Ar_orth*y_1;.. P_Ar_orth*y_M]||^2/(MN-MQ)
        for q = 1:param.numInt
            b_q = param.B_RS(:,q); % the q-th column of B
            at_int_est_n_l_q = reshape((param.Ar_int*b_q)'*y_rD_3D_n_l,[],1); % estimator of at_int_q = [(param.Ar_int*b_q)'*y_1; ...; (param.Ar_int*b_q)'*y_M]
            b_tuta_estRaw(n,l,:,q) = fft(at_int_est_n_l_q,param.Nt)/param.Nt;
        end
        hsquare_estRaw(n,l,:,:) = abs(squeeze(b_tuta_estRaw(n,l,:,:))).^2;
        einr_estRaw(n,l,:,:) = hsquare_estRaw(n,l,:,:)/noiseEstRaw(n,l);
    end
end
% EINR estimation
for n = 1: param.Nrfft
    for l = 1: param.Nvfft
        einr_est(n,l,:,:) = getEinrEstSingleBin(n,l,hsquare_estRaw,noiseEstRaw,param);
    end
end
end
