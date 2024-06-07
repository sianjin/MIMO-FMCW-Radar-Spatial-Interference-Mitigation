function [at_int_estClairvoyant,param] = getEstIntTxSteeringVec(X_3D_Int,param)
%% Range FFT
% Reference: Basic spectral analysis
% Web: https://www.mathworks.com/help/matlab/math/basic-spectral-analysis.html
% Range FFT of the combined signals from multiple Txs
% As the phases of the received signals from Txs are different,
% they might add constructively or destructively
Y_r_3D  = complex(zeros(param.Nrfft,param.Nr,param.Nsweep));
for l = 1: param.Nsweep
    for m = 1:param.Nr
        X_3D_m_l = X_3D_Int(:,m,l);
        Y_r_3D(:,m,l) = fftshift(fft(X_3D_m_l,param.Nrfft));
    end
end
%% Co-located Tx phase decoding
Y_vm_r_3D_raw = complex(zeros([size(Y_r_3D),param.Nt]));
% Phase decoding
for l = 1: param.Nsweep
    % Phase decoding on each element
    for mt = 1: param.Nt
        wl_mt = param.w(mt,l); % No need to take conjugate because dechirping already conjugates target signal
        Y_vm_r_3D_raw(:,:,l,mt) = Y_r_3D(:,:,l)*wl_mt;
    end
end
Y_v_r_3D = complex(zeros(param.Nrfft,param.Nt*param.Nr,param.Nsweep));
if param.slowTimeCodeType == "TDM"
    param.Nvfft = 2^nextpow2(param.Nsweep); % length of velocity FFT
    param.Tg = param.Tg*param.Nt;
    param.vmax = param.vmax/param.Nt;
    param.binVelococityTarget = ceil(param.vTarget/(param.c/param.fc/2)/(-1/param.Tg)*param.Nvfft + param.Nvfft/2); % velocity bin of the target
    Y_vm_r_3D_rawNew = complex(zeros(param.Nrfft,param.Nr,param.Nsweep/param.Nt,param.Nt));
    for n = 1: param.Nrfft
        for mr = 1:param.Nr
            for mt = 1:param.Nt
                Y_vm_r_3D_raw_n_mr_mt = Y_vm_r_3D_raw(n,mr,:,mt);
                Y_vm_r_3D_raw_n_mr_mt(Y_vm_r_3D_raw_n_mr_mt==0) = [];
                Y_vm_r_3D_rawNew(n,mr,:,mt) = Y_vm_r_3D_raw_n_mr_mt;
            end
        end
    end
    Y_vm_r_3D_raw = Y_vm_r_3D_rawNew;
    Y_v_r_3D = complex(zeros(param.Nrfft,param.Nt*param.Nr,param.Nsweep/param.Nt));
end
for mt = 1:param.Nt
    Y_v_r_3D(:,((mt-1)*param.Nr+1):(mt*param.Nr),:) = Y_vm_r_3D_raw(:,:,:,mt);
end
%% Veclocity FFT for phase coded MIMO radar
param.Nvfft = 2^nextpow2(param.Nsweep); % length of velocity FFT
Y_rD_3D = complex(zeros(param.Nrfft,param.Nr*param.Nt,param.Nvfft));
% Y_vm_rD_3D = complex(zeros(param.Nrfft,param.Nr,param.Nvfft,param.Nt));
for n = 1: param.Nrfft
    for m = 1:param.Nr*param.Nt
        Y_r_3D_n_m = Y_v_r_3D(n,m,:); 
        Y_rD_3D(n,m,:) = fftshift(fft(Y_r_3D_n_m,param.Nvfft));
    end
end
%% Angle FFT for phase coded MIMO radar
param.Nafft = 2^nextpow2(param.Nr*param.Nt); % length of angle FFT
param.Nacorr = param.Nafft;
% Interference Rx steering vector
int_freq = sind(param.azIntVec)/param.lambda*param.rxEleSpacing;
param.Ar_int = exp(1j*2*pi*(0:param.Nr-1)'*int_freq);
B_RS_preinv = param.Ar_int'*param.Ar_int;
param.B_RS = inv(B_RS_preinv); 
param.P_Ar_int_orth = eye(param.Nr) - param.Ar_int*(B_RS_preinv\param.Ar_int');
at_int_estClairvoyant = zeros(param.Nt,param.numInt); % esimation of at_int on different range-velocity bins 
Y_rD_3D_n_l = Y_rD_3D(param.binRangeTarget(param.targetIdx),:,param.binVelococityTarget(param.targetIdx))';
% Obtain estimator of decoded interference Tx steering vector at_int
y_rD_3D_n_l = reshape(Y_rD_3D_n_l,[param.Nr,param.Nt]); % decomponse y into y_1, y_2, ..., y_M
for q = 1:param.numInt
    b_q = param.B_RS(:,q); % the q-th column of B
    at_int_est_q = reshape((param.Ar_int*b_q)'*y_rD_3D_n_l,[],1); % estimator of at_int_q = [(param.Ar_int*b_q)'*y_1; ...; (param.Ar_int*b_q)'*y_M]
    at_int_estClairvoyant(:,q) = at_int_est_q;
end
% Validation of estimator of at_int_q: Y_rD_3D_n_l approximately = sum_q kron(at_int_est_q,param.Ar_int_q)
% Y_rD_3D_n_l_est = kron(at_int_est_record(n,:,l,1).',param.Ar_int(:,1)) + kron(at_int_est_record(n,:,l,2).',param.Ar_int(:,2));
end
