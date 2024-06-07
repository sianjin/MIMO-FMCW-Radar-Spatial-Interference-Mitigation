% Obtain EINR estimation for each range-velocity bins
function R_est = getIntCovEst(Y_rD_3D,param)
noiseEstRaw = zeros(param.Nrfft,param.Nvfft); % esimation of noise power
at_int_est = zeros(param.Nrfft,param.Nvfft,param.Nt,param.numInt);
R_est = zeros(param.Nrfft,param.Nvfft,param.Nt*param.Nr,param.Nt*param.Nr); % esimation of EINR from nearby range-velocity bins 
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
            at_int_est(n,l,:,q) = reshape((param.Ar_int*b_q)'*y_rD_3D_n_l,[],1); % estimator of at_int_q = [(param.Ar_int*b_q)'*y_1; ...; (param.Ar_int*b_q)'*y_M]
        end
    end
end
% EINR estimation
for n = 1: param.Nrfft
    for l = 1: param.Nvfft
        R_est(n,l,:,:) = getIntCovEstSingleBin(n,l,at_int_est,noiseEstRaw,param);
    end
end
end
