function param = getDetectResult(param)
param.Nafft = 2^nextpow2(param.Nr*param.Nt); % length of angle FFT
param.Nrfft = 2^(nextpow2(ceil(param.FsADC*param.Tg))-1); % length of range FFT
param.fftIntPowSet = zeros(1,param.numIterSim);
param.fastClairvoyantIntPowSet = zeros(1,param.numIterSim);
param.fastRsIntPowSet = zeros(1,param.numIterSim);
param.fastGsIntPowSet = zeros(1,param.numIterSim);
param.fastLcmvIntPowSet = zeros(1,param.numIterSim);
param.agsIntPowSet = zeros(1,param.numIterSim);
for iterSim = 1: param.numIterSim
%% Dechirped mixed target and interference signal
param.dechirpMix = param.dechirpTar;
param.azIntVec = -80+rand*160; % Interference angle
for numIntIdx = 1:param.numInt
    % Setup unique parameters for interference]
    param.dInt = ceil((param.MinIntDist+2*rand)*10)/10; % interference distance, m
    param.vInt = ceil((-param.vmax+2*param.vmax*rand)*10)/10; % real int velocity vInt, detected int velocity=1/2*vInt
    param.azInt = param.azIntVec(numIntIdx); % interference azimuth angle, degree 
    param.intMotion = phased.Platform('InitialPosition',[param.dInt.*cosd(param.azInt);param.dInt.*sind(param.azInt);0],...
        'Velocity',[-param.vInt;0;0]);
    param.hcInt = (10+10*rand)*1e12; % interference chirp sloped
    param.TcInt = ceil(param.Bc/param.hcInt*1e7)/1e7; % interference chirp duration
    param.TiInt = (6+2*rand)/1e6; % idle inter-chirp duration
    param.TgInt = param.TcInt + param.TiInt; % interference chirp repetition period
    param.toffInt = param.Tc*rand; % initial time offset between interference 1 and victim radar
    dechirpInt = getDechirpedInt(param);
    param.dechirpMix = param.dechirpMix + dechirpInt;
end
%% Get lowpassed and ADC sampled data
[X_3D,param] = getRawADCData(param);
[X_3D_Int,param] = getRawADCDataInt(param);
%% RAV Processing setup
param.Nvfft = 2^nextpow2(param.Nsweep); % length of velocity FFT
param.Nafft = 2^nextpow2(param.Nr*param.Nt); % length of angle FFT
param.binRangeTarget = ceil(param.dTarget/(param.c/2/param.hc)/param.FsADC*param.Nrfft + param.Nrfft/2+1)+3; % range bin of the target
param.binVelococityTarget = ceil(param.vTarget/(param.c/param.fc/2)/(-1/param.Tg)*param.Nvfft + param.Nvfft/2)+1; % velocity bin of the target
param.binAngleTarget = ceil(sind(param.azTarget)/param.lambda*param.rxEleSpacing*param.Nafft + param.Nafft/2);  % angle bin of the target
param.binAngleInt = ceil(sind(param.azIntVec)/param.lambda*param.rxEleSpacing*param.Nafft + param.Nafft/2);  % angle bin of the interference
%% Get interference Tx steering vector
[at_int_estClairvoyant,param] = getEstIntTxSteeringVec(X_3D_Int,param);
%% Range FFT
% Reference: Basic spectral analysis
% Web: https://www.mathworks.com/help/matlab/math/basic-spectral-analysis.html
% Range FFT of the combined signals from multiple Txs
% As the phases of the received signals from Txs are different,
% they might add constructively or destructively
Y_r_3D  = complex(zeros(param.Nrfft,param.Nr,param.Nsweep));
for l = 1: param.Nsweep
    for m = 1:param.Nr
        X_3D_m_l = X_3D(:,m,l);
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
if param.slowTimeCodeType == "TDM"
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
    param.Nsweep = param.Nsweep/param.Nt;
    Y_vm_r_3D_raw = Y_vm_r_3D_rawNew;
end
Y_v_r_3D = complex(zeros(param.Nrfft,param.Nt*param.Nr,param.Nsweep));
for mt = 1:param.Nt
    Y_v_r_3D(:,((mt-1)*param.Nr+1):(mt*param.Nr),:) = Y_vm_r_3D_raw(:,:,:,mt);
end
%% Veclocity FFT for phase coded MIMO radar
Y_rD_3D = complex(zeros(param.Nrfft,param.Nr*param.Nt,param.Nvfft));
% Y_vm_rD_3D = complex(zeros(param.Nrfft,param.Nr,param.Nvfft,param.Nt));
for n = 1: param.Nrfft
    for m = 1:param.Nr*param.Nt
        Y_r_3D_n_m = Y_v_r_3D(n,m,:); 
        Y_rD_3D(n,m,:) = fftshift(fft(Y_r_3D_n_m,param.Nvfft));
    end
end
%% Angle detection for phase coded MIMO radar
fftStats = complex(zeros(param.Nrfft,param.Nafft));
fastClairvoyantStats = complex(zeros(param.Nrfft,param.Nafft));
fastRsStats = complex(zeros(param.Nrfft,param.Nafft));
fastGsStats = complex(zeros(param.Nrfft,param.Nafft));
fastLcmvStats = complex(zeros(param.Nrfft,param.Nafft));
agsStats = complex(zeros(param.Nrfft,param.Nafft));
intStatRangeBin =param.binRangeTarget(param.targetIdx);
intStatVelocityBin = param.binVelococityTarget(param.targetIdx);
intStatAngleBin = param.binAngleTarget(param.targetIdx);
% Noise and covariance estimation parameters
rangeGI = 8; % gurad interval on range domain
velocityGI = 4; % gurad interval on velocity domain
rangeTrainCells = 4; % number of training cells on each side of range domain
velocityTrainCells = 4; % number of training cells on each side of velocity domain
% Fast RS: calculate denominator
int_freq = sind(param.azIntVec)/param.lambda*param.rxEleSpacing;
param.Ar_int = exp(1j*2*pi*(0:param.Nr-1)'*int_freq);
B_RS_preinv = param.Ar_int'*param.Ar_int;
param.P_Ar_int_orth = eye(param.Nr) - param.Ar_int*(B_RS_preinv\param.Ar_int');
param.B_RS = inv(B_RS_preinv); 
[U_RS,D_RS] = eig(param.B_RS); % eigenvalue decomposition of (Ar_int'*Ar_int)^{-1}
d_RS = abs(diag(D_RS)); % a column vector of eigenvalues of (Ar_int'*Ar_int)^{-1}
U_tuta_RS = param.Ar_int*U_RS;
den_rs = param.Nr;
for q = 1:param.numInt
    den_rs = den_rs - d_RS(q)*abs(fftshift(fft([U_tuta_RS(:,q);zeros(param.Nr*param.Nt-param.Nr,1)],param.Nr*param.Nt))).^2;
end
den_rs = den_rs*param.Nt;
% Detection result on each range bin
for n = 1: param.Nrfft
    Y_rD_3D_n = Y_rD_3D(n,:,intStatVelocityBin)';
    % Angle FFT
    fft_y_n = fftshift(fft(Y_rD_3D_n,param.Nafft));
    fftStats(n,:) = abs(fft_y_n).^2/(param.Nt*param.Nr);
    % Fast clairvoyant
    b_tuta_n = zeros(param.Nt,param.numInt);
    for q = 1:param.numInt
        at_int_n_q = at_int_estClairvoyant(n,:,q);
        b_tuta_n(:,q) = fft(at_int_n_q,param.Nt)/param.Nt;
    end
    clairvoyant_n = fft_y_n;
    for q = 1:param.numInt
        ar_int_q = param.Ar_int(:,q); % q-th interference Rx steering vector
        clairvoyant_int_q = param.Nt*kron(ones(param.Nr,1),b_tuta_n(:,q)).*fftshift(fft([ar_int_q;zeros(param.Nr*param.Nt-param.Nr,1)],param.Nr*param.Nt));
        clairvoyant_n = clairvoyant_n - clairvoyant_int_q;
    end
    fastClairvoyantStats(n,:) = abs(clairvoyant_n).^2/(param.Nt*param.Nr);
    % Fast RS
    y_rD_3D_n = reshape(Y_rD_3D_n,[param.Nr,param.Nt]); % decomponse y into y_1, y_2, ..., y_M
    Y_rD_3D_n_proj = reshape(param.P_Ar_int_orth*y_rD_3D_n,[],1);
    rs_n = fftshift(fft(Y_rD_3D_n_proj,param.Nafft));
    num_rs = abs(rs_n).^2;
    fastRsStats(n,:) = num_rs./den_rs;
    % Fast GS
    fast_einr_est = getEinrEstFast(Y_rD_3D,n,intStatVelocityBin,param,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells);
    Lambda_einr = zeros(param.numInt,param.numInt,param.Nt);
    for mt = 1:param.Nt
        Lambda_einr(:,:,mt) = diag(fast_einr_est(mt,:));
    end
    den_gs = zeros(param.Nr,param.Nt);
    P_Ar_int_orth_reg = cell(param.Nt,1);
    for mt = 1:param.Nt
        B_GS_preinv_mt = inv(Lambda_einr(:,:,mt))/param.Nt + param.Ar_int'*param.Ar_int;
        B_GS_mt = inv(B_GS_preinv_mt);
        [U_GS_mt,D_GS_mt] = eig(B_GS_mt);
        U_tuta_GS_mt = param.Ar_int*U_GS_mt;
        P_Ar_int_orth_reg{mt} = eye(param.Nr) - U_tuta_GS_mt*D_GS_mt*U_tuta_GS_mt';
        d_GS_mt = abs(diag(D_GS_mt));
        den_gs_mt = param.Nr;
        rotaVec_Nr_mt = exp(-1j*2*pi*(mt-1)*(0:param.Nr-1)'/(param.Nr*param.Nt));
        for q = 1:param.numInt
            den_gs_mt = den_gs_mt - ...
                d_GS_mt(q)*abs(fftshift(fft([U_tuta_GS_mt(:,q)].*rotaVec_Nr_mt,param.Nr))).^2;
        end
        den_gs(:,mt) = den_gs_mt*param.Nt;
    end
    den_gs = reshape(den_gs',[],1);
    num_gs = zeros(param.Nr,param.Nt);
    for mt = 1:param.Nt
        P_Ar_int_orth_reg_mt = P_Ar_int_orth_reg{mt};
        Y_rD_3D_n_proj_mt = P_Ar_int_orth_reg_mt*y_rD_3D_n;
        Y_rD_3D_n_proj_mt_vec = reshape(Y_rD_3D_n_proj_mt,[],1);
        num_gs_full_mt = fftshift(fft(Y_rD_3D_n_proj_mt_vec,param.Nafft));
        num_gs(:,mt) = num_gs_full_mt(mt:param.Nt:end);
    end
    num_gs = abs(num_gs).^2;
    num_gs = reshape(num_gs',[],1);
    fastGsStats(n,:) = num_gs./den_gs;
    % Fast LCMV
    fast_R_est = getIntCovEstFast(Y_rD_3D,n,intStatVelocityBin,param,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells);
    R_est_inv = inv(fast_R_est);
    [U_LCMV,D_LCMV] = eig(R_est_inv); % eigenvalue decomposition of (R_est_n_l)^{-1}
    d_LCMV = abs(diag(D_LCMV)); % a column vector of eigenvalues of (R_est_n_l)^{-1}
    den_lcmv = zeros(param.Nr*param.Nt,1);
    for m = 1:param.Nr*param.Nt
        den_lcmv = den_lcmv + d_LCMV(m)*abs(fftshift(fft(U_LCMV(:,m),param.Nr*param.Nt))).^2;
    end
    Y_rD_3D_n_whiten = R_est_inv*Y_rD_3D_n;
    lcmv_n = fftshift(fft(Y_rD_3D_n_whiten,param.Nafft));
    num_lcmv = abs(lcmv_n).^2;
    fastLcmvStats(n,:) = num_lcmv./den_lcmv;
    % AGS
    for m = [intStatAngleBin,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft)]
        at = exp(1j*2*pi*(0:param.Nt-1)'*param.txEleSpacing/param.lambda*(-1+(m-1)*2/param.Nafft));
        ar = exp(1j*2*pi*(0:param.Nr-1)'*param.rxEleSpacing/param.lambda*(-1+(m-1)*2/param.Nafft));
        av = kron(at,ar);
        R_recon = getCovRecon(fast_R_est,param.Nt,param.Nr,at);
        agsStats(n,m) = abs((R_recon\av)'*Y_rD_3D_n)^2/abs(av'*(R_recon\av));
    end
end
%% Plot
% param.PowAngleFFTdB = 10*log10(fftStats);
% param.PowClairvoyantdetectStatsdB = 10*log10(fastClairvoyantStats);
% tarClairvoyantdetectStatsdB = param.PowClairvoyantdetectStatsdB(intStatRangeBin,intStatAngleBin);
% param.PowRSStatsdB = 10*log10(fastRsStats);
% param.PowGSStatsdB = 10*log10(fastGsStats);
% param.PowLCMVStatsdB = 10*log10(fastLcmvStats) + tarClairvoyantdetectStatsdB - 10*log10(fastLcmvStats(intStatRangeBin,intStatAngleBin));
% param.PowAGSStatsdB = 10*log10(agsStats) + tarClairvoyantdetectStatsdB - 10*log10(agsStats(intStatRangeBin,intStatAngleBin));
% %% Plot RA map for angle FFT
% figure
% f_r = param.FsADC*(-param.Nrfft/2:param.Nrfft/2-1)/param.Nrfft; % Hz
% range = f_r * param.c/2/param.hc;
% tarAngle = asind(param.lambda/param.rxEleSpacing*(-(param.Nafft/2):(param.Nafft/2)-1)/param.Nafft); % Hz
% [Range, Angle] = meshgrid(range, tarAngle);
% mesh(Range, Angle, param.PowAngleFFTdB.') 
% colorbar
% xlim([0, param.dmax])
% xlabel('Range (m)')
% ylabel('Angle (degree)')
% zlabel('Power (dB)')
% title('Range-Angle Image, Angle FFT')
% %% Plot RA map for clairvoyant
% figure
% [Range, Angle] = meshgrid(range, tarAngle);
% mesh(Range, Angle, param.PowClairvoyantdetectStatsdB.') 
% colorbar
% xlim([0, param.dmax])
% xlabel('Range (m)')
% ylabel('Angle (degree)')
% zlabel('Power (dB)')
% title('Range-Angle Image, Clairvoyant')
% %% Plot RA map for RS
% figure
% [Range, Angle] = meshgrid(range, tarAngle);
% mesh(Range, Angle, param.PowRSStatsdB.') 
% colorbar
% xlim([0, param.dmax])
% xlabel('Range (m)')
% ylabel('Angle (degree)')
% zlabel('Power (dB)')
% title('Range-Angle Image, RS')
% %% Plot RA map for GS
% figure
% [Range, Angle] = meshgrid(range, tarAngle);
% mesh(Range, Angle, param.PowGSStatsdB.')
% colorbar
% xlim([0, param.dmax])
% xlabel('Range (m)')
% ylabel('Angle (degree)')
% zlabel('Power (dB)')
% title('Range-Angle Image, GS')
% %% Plot RA map for LCMV-SMI
% figure
% [Range, Angle] = meshgrid(range, tarAngle);
% mesh(Range, Angle, param.PowLCMVStatsdB.')
% colorbar
% xlim([0, param.dmax])
% xlabel('Range (m)')
% ylabel('Angle (degree)')
% zlabel('Power (dB)')
% title('Range-Angle Image, LCMV-SMI')
% %% Plot RA map for AGS
% figure
% [Range, Angle] = meshgrid(range, tarAngle);
% mesh(Range, Angle, param.PowAGSStatsdB.')
% colorbar
% xlim([0, param.dmax])
% xlabel('Range (m)')
% ylabel('Angle (degree)')
% zlabel('Power (dB)')
% title('Range-Angle Image, AGS')
%% RA-domain SINR
% Extended number of target range bins on both sides 
binRangeExtend = 14;
% Extended target angle bins on both sides  
binAnlgeExtend = 1;
% Angle FFT
fftNonTargetAmplitude = fftStats;
for tarIdx = 1:length(param.dTarget)
    binRangeTargetExtendTarIdx = max(param.binRangeTarget(tarIdx)-binRangeExtend,1):min(param.binRangeTarget(tarIdx)+binRangeExtend,param.Nrfft); 
    binAngleTargetExtendTarIdx = max(param.binAngleTarget(tarIdx)-binAnlgeExtend,1)-1:min(param.binAngleTarget(tarIdx)+binAnlgeExtend,param.Nafft); 
    fftNonTargetAmplitude(binRangeTargetExtendTarIdx,binAngleTargetExtendTarIdx) = 0;
end
fftNonTargetPw = fftNonTargetAmplitude.^2;
fftIntPw = fftNonTargetPw(:,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft));
fftIntPw = fftIntPw(fftIntPw~=0);
fftIntAvgPw = mean(fftIntPw,'all');
param.fftIntPowSet(iterSim) = pow2db(fftIntAvgPw);
% Clairvoyant
clairvoyantNonTargetAmplitude = fastClairvoyantStats;
for tarIdx = 1:length(param.dTarget)
    binRangeTargetExtendTarIdx = max(param.binRangeTarget(tarIdx)-binRangeExtend,1):min(param.binRangeTarget(tarIdx)+binRangeExtend,param.Nrfft); 
    binAngleTargetExtendTarIdx = max(param.binAngleTarget(tarIdx)-binAnlgeExtend,1)-1:min(param.binAngleTarget(tarIdx)+binAnlgeExtend,param.Nafft); 
    clairvoyantNonTargetAmplitude(binRangeTargetExtendTarIdx,binAngleTargetExtendTarIdx) = 0;
end
clairvoyantNonTargetPw = clairvoyantNonTargetAmplitude.^2;
clairvoyantIntPw = clairvoyantNonTargetPw(:,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft));
clairvoyantIntPw = clairvoyantIntPw(clairvoyantIntPw~=0);
clairvoyantIntAvgPw = mean(clairvoyantIntPw,'all');
param.fastClairvoyantIntPowSet(iterSim) = pow2db(clairvoyantIntAvgPw);
% RS
rsNonTargetAmplitude = fastRsStats;
for tarIdx = 1:length(param.dTarget)
    binRangeTargetExtendTarIdx = max(param.binRangeTarget(tarIdx)-binRangeExtend,1):min(param.binRangeTarget(tarIdx)+binRangeExtend,param.Nrfft); 
    binAngleTargetExtendTarIdx = max(param.binAngleTarget(tarIdx)-binAnlgeExtend,1)-1:min(param.binAngleTarget(tarIdx)+binAnlgeExtend,param.Nafft); 
    rsNonTargetAmplitude(binRangeTargetExtendTarIdx,binAngleTargetExtendTarIdx) = 0;
end
rsNonTargetPw = rsNonTargetAmplitude.^2;
rsIntPw = rsNonTargetPw(:,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft));
rsIntPw = rsIntPw(rsIntPw~=0);
rsIntAvgPw = mean(rsIntPw,'all');
param.fastRsIntPowSet(iterSim) = pow2db(rsIntAvgPw);
% GS
gsNonTargetAmplitude = fastGsStats;
for tarIdx = 1:length(param.dTarget)
    binRangeTargetExtendTarIdx = max(param.binRangeTarget(tarIdx)-binRangeExtend,1):min(param.binRangeTarget(tarIdx)+binRangeExtend,param.Nrfft); 
    binAngleTargetExtendTarIdx = max(param.binAngleTarget(tarIdx)-binAnlgeExtend,1)-1:min(param.binAngleTarget(tarIdx)+binAnlgeExtend,param.Nafft); 
    gsNonTargetAmplitude(binRangeTargetExtendTarIdx,binAngleTargetExtendTarIdx) = 0;
end
gsNonTargetPw = gsNonTargetAmplitude.^2;
gsIntPw = gsNonTargetPw(:,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft));
gsIntPw = gsIntPw(gsIntPw~=0);
gsIntAvgPw = mean(gsIntPw,'all');
param.fastGsIntPowSet(iterSim) = pow2db(gsIntAvgPw);
% LCMV-SMI
lcmvNonTargetAmplitude = fastLcmvStats/fastLcmvStats(intStatRangeBin,intStatAngleBin)*fastClairvoyantStats(intStatRangeBin,intStatAngleBin);
for tarIdx = 1:length(param.dTarget)
    binRangeTargetExtendTarIdx = max(param.binRangeTarget(tarIdx)-binRangeExtend,1):min(param.binRangeTarget(tarIdx)+binRangeExtend,param.Nrfft); 
    binAngleTargetExtendTarIdx = max(param.binAngleTarget(tarIdx)-binAnlgeExtend,1)-1:min(param.binAngleTarget(tarIdx)+binAnlgeExtend,param.Nafft); 
    lcmvNonTargetAmplitude(binRangeTargetExtendTarIdx,binAngleTargetExtendTarIdx) = 0;
end
lcmvNonTargetPw = lcmvNonTargetAmplitude.^2;
lcmvIntPw = lcmvNonTargetPw(:,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft));
lcmvIntPw = lcmvIntPw(lcmvIntPw~=0);
lcmvIntAvgPw = mean(lcmvIntPw,'all');
param.fastLcmvIntPowSet(iterSim) = pow2db(lcmvIntAvgPw);
% AGS
agsNonTargetAmplitude = agsStats/agsStats(intStatRangeBin,intStatAngleBin)*fastClairvoyantStats(intStatRangeBin,intStatAngleBin);
for tarIdx = 1:length(param.dTarget)
    binRangeTargetExtendTarIdx = max(param.binRangeTarget(tarIdx)-binRangeExtend,1):min(param.binRangeTarget(tarIdx)+binRangeExtend,param.Nrfft); 
    binAngleTargetExtendTarIdx = max(param.binAngleTarget(tarIdx)-binAnlgeExtend,1)-1:min(param.binAngleTarget(tarIdx)+binAnlgeExtend,param.Nafft); 
    agsNonTargetAmplitude(binRangeTargetExtendTarIdx,binAngleTargetExtendTarIdx) = 0;
end
agsNonTargetPw = agsNonTargetAmplitude.^2;
agsIntPw = agsNonTargetPw(:,max(param.binAngleInt-2,1):min(param.binAngleInt+2,param.Nafft));
agsIntPw = agsIntPw(agsIntPw~=0);
agsIntAvgPw = mean(agsIntPw,'all');
param.agsIntPowSet(iterSim) = pow2db(agsIntAvgPw);
end
end