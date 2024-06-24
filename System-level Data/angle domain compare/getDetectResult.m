function param = getDetectResult(param)
%% Dechirped mixed target and interference signal
param.dechirpMix = param.dechirpTar;
for numIntIdx = 1:param.numInt
    % Setup unique parameters for interference]
    param.dInt = ceil((param.AvgDist+2*rand)*10)/10; % interference distance, m
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
%% RV Processing setup
param.Nvfft = 2^nextpow2(param.Nsweep); % length of velocity FFT
param.binRangeTarget = ceil(param.dTarget/(param.c/2/param.hc)/param.FsADC*param.Nrfft + param.Nrfft/2+1)+1; % range bin of the target
param.binVelococityTarget = ceil(param.vTarget/(param.c/param.fc/2)/(-1/param.Tg)*param.Nvfft + param.Nvfft/2); % velocity bin of the target
%% Get interference Tx steering vector
[at_int_est_recordClairvoyant,param] = getEstIntTxSteeringVec(X_3D_Int,param);
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
param.Nvfft = 2^nextpow2(param.Nsweep); % length of velocity FFT
Y_rD_3D = complex(zeros(param.Nrfft,param.Nr*param.Nt,param.Nvfft));
% Y_vm_rD_3D = complex(zeros(param.Nrfft,param.Nr,param.Nvfft,param.Nt));
for n = 1: param.Nrfft
    for m = 1:param.Nr*param.Nt
        Y_r_3D_n_m = Y_v_r_3D(n,m,:); 
        Y_rD_3D(n,m,:) = fftshift(fft(Y_r_3D_n_m,param.Nvfft));
    end
end
%% Plot RV map for BPM MIMO radar
% figure
% f_r = param.FsADC*(-param.Nrfft/2:param.Nrfft/2-1)/param.Nrfft; % Hz
% range = f_r * param.c/2/param.hc;
% f_D = -1/(param.Tg)*(-(param.Nvfft/2):(param.Nvfft/2)-1)/param.Nvfft; % Hz
% velocity = f_D * param.c/param.fc/2;
% [Range, Velocity] = meshgrid(range, velocity);
% % Non-coherent integration
% Amp_Y_rD_3D_linear = abs(Y_rD_3D);
% Pow_Y_rD_3D_linear = Amp_Y_rD_3D_linear.^2;
% mesh(Range, Velocity, 10*log10(squeeze(sum(Pow_Y_rD_3D_linear,2)).')) % 2D-FFT result
% colorbar
% colormap(hot)
% cmap = colormap;
% cmapNew = flip(cmap);
% colormap(cmapNew)
% xlim([0, param.dmax])
% % As we adopt pseudo random code MIMO without filtering, the max velocity
% % is unchanged
% ylim([-param.vmax,param.vmax])
% xlabel('Range (m)')
% ylabel('Velocity (m/s)')
% zlabel('Power (dB)')
% title('Range-Velocity Image')
%% RD CFAR Detection
% countY_rD_2D = zeros(param.Nrfft,param.Nvfft);
% detector = phased.CFARDetector2D('TrainingBandSize',[5,5], ...
%     'ThresholdFactor','Auto','GuardBandSize',[3,3], ...
%     'ProbabilityFalseAlarm',1e-15,'ThresholdOutputPort',true);
% Ngc = detector.GuardBandSize(2);
% Ngr = detector.GuardBandSize(1);
% Ntc = detector.TrainingBandSize(2);
% Ntr = detector.TrainingBandSize(1);
% cutidx = [];
% colstart = Ntc + Ngc + 1;
% colend = param.Nvfft - ( Ntc + Ngc);
% rowstart = Ntr + Ngr + 1;
% rowend = param.Nrfft - ( Ntr + Ngr);
% for m = colstart:colend
%     for n = rowstart:rowend
%         cutidx = [cutidx,[n;m]];
%     end
% end
% dets = detector(squeeze(sum(Pow_Y_rD_3D_linear,2)),cutidx);
% for k = 1 : length(dets)
%     countY_rD_2D(cutidx(1,k), cutidx(2,k)) = dets(k);
% end
%% Plot RV detection result for phase coded MIMO radar
% figure
% mesh(Range, Velocity, countY_rD_2D.') % cfar detection result
% colorbar
% colormap(hot)
% cmap = colormap;
% cmapNew = flip(cmap);
% colormap(cmapNew)
% xlim([0, param.dmax])
% % As we adopt pseudo random code MIMO without filtering, the max velocity
% % is unchanged
% ylim([-param.vmax,param.vmax])
% zlim([0, 1]);
% xlabel('Range (m)')
% ylabel('Velocity (m/s)')
% title('Detection Result on Range-Velocity Bins')
%% Angle spectrum for phase coded MIMO radar
param.Nafft = 2^nextpow2(param.Nr*param.Nt*4); % length of angle FFT
clairvoyantStats = complex(zeros(param.Nafft,1));
lcmvStats = complex(zeros(param.Nafft,1));
agsStats = complex(zeros(param.Nafft,1));
% Range-velocity bins for fast interference statistic estimation
intStatRangeBin = param.binRangeTarget(param.targetIdx);
intStatVelocityBin = param.binVelococityTarget(param.targetIdx);
% Noise and covariance estimation parameters
rangeGI = 8; % range domain gurad interval on each side
velocityGI = 4; % velocity domain gurad interval on each side
rangeTrainCells = 4; % number of training cells on each side of range domain
velocityTrainCells = 4; % number of training cells on each side of velocity domain
% Get interference covariance matrix estimation from nearby range-velocity
% bins for LCMV detector
% R_est = getIntCovEst(Y_rD_3D,param);
fast_R_est = getIntCovEstFast(Y_rD_3D,intStatRangeBin,intStatVelocityBin,param,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells);
% Obtain virtual array signal y on range bin n and velocity bin l
Y_rD_3D_n_l = Y_rD_3D(intStatRangeBin,:,intStatVelocityBin)';
Y_rD_3D_n_l = [Y_rD_3D_n_l;zeros(param.Nafft-param.Nr*param.Nt,1)]; % zero-pad virtual signal
% Obtain estimator of decoded interference Tx steering vector
% at_int_q, q=1,2,...,Q, and essential interference amplitude
% b_tutua_q, q=1,2,...,Q, at range bin n and velocity bin l
At_int_n_l = zeros(param.Nt,param.numInt);
for q = 1:param.numInt
    At_int_n_l(:,q) = at_int_est_recordClairvoyant(:,q);
end
% Obtain EINR
% einr = squeeze(einr_est(n,l,:,:));
einr = fast_einr_est;
Lambda_einr = zeros(param.numInt,param.numInt,param.Nt);
for mt = 1:param.Nt
    Lambda_einr(:,:,mt) = diag(einr(mt,:));
end
% Angle FFT
Y_raD_3D_n_l_fft = fftshift(fft(Y_rD_3D_n_l,param.Nafft));
fftStats = abs(Y_raD_3D_n_l_fft).^2/(param.Nt*param.Nr);
for m = 1:param.Nafft
    % Steering vector swept over each angle bin
    at = exp(1j*2*pi*(0:param.Nt-1)'*param.txEleSpacing/param.lambda*(-1+(m-1)*2/param.Nafft));
    ar = exp(1j*2*pi*(0:param.Nr-1)'*param.rxEleSpacing/param.lambda*(-1+(m-1)*2/param.Nafft));
    av = kron(at,ar);
    % Clairvoyant
    av_zeroPad_clairvoyant = [av;zeros(param.Nafft-param.Nr*param.Nt,1)];
    clairvoyantStats_n_m_l = av_zeroPad_clairvoyant'*Y_rD_3D_n_l;
    for q = 1:param.numInt
        ar_int_q = param.Ar_int(:,q); % q-th interference Rx steering vector
        b_tuta_n_l_q = at'*At_int_n_l(:,q)/norm(at)^2;
        clairvoyantStats_n_m_l = clairvoyantStats_n_m_l-b_tuta_n_l_q*param.Nt*ar'*ar_int_q;
    end
    clairvoyantStats(m) = abs(clairvoyantStats_n_m_l)^2/(param.Nt*param.Nr);
    % LCMV
    R_est_n_l = fast_R_est;
    av_zeroPad_lcmv = [(R_est_n_l\av);zeros(param.Nafft-param.Nr*param.Nt,1)];
    lcmvStats(m) = abs(av_zeroPad_lcmv'*Y_rD_3D_n_l)^2/abs(av'*(R_est_n_l\av));
    % AGS
    R_recon_n_l = getCovRecon(R_est_n_l,param.Nt,param.Nr,at);
    av_zeroPad_ags = [(R_recon_n_l\av);zeros(param.Nafft-param.Nr*param.Nt,1)];
    agsStats(m) = abs(av_zeroPad_ags'*Y_rD_3D_n_l)^2/abs(av'*(R_recon_n_l\av));
end
param.PowAngleFFTdB = 10*log10(fftStats);
param.PowClairvoyantdetectStatsdB = 10*log10(clairvoyantStats);
[maxPowClairvoyantdetectStatsdB,maxClairIdx] = max(param.PowClairvoyantdetectStatsdB);
param.PowLCMVStatsdB = 10*log10(lcmvStats) + maxPowClairvoyantdetectStatsdB - 10*log10(lcmvStats(maxClairIdx));
param.PowAGSStatsdB = 10*log10(agsStats) + maxPowClairvoyantdetectStatsdB - 10*log10(agsStats(maxClairIdx));
