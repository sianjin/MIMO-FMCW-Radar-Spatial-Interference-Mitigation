function [X_3D,param] = getRawADCData(param)
%% Low pass filter the IF signal
X_3D_fullsample = complex(zeros(param.Ng,param.Nr,param.Nsweep)); % dechirped and lowpassed data cube
[~,digFilter] = lowpass(param.dechirpMix(:,1,1),param.fH,param.Fs,'Steepness',0.999,...
    'ImpulseResponse','fir','StopbandAttenuation',80); % pass received signal into LPF
N = filtord(digFilter);
% fvtool(digFilter,'Analysis','freq')
for l = 1:param.Nsweep
    for m = 1:param.Nr
        dechirpMixl = param.dechirpMix(:,m,l); % dechirp received target signal
        % Compensate for the delay introduced by the filter
        filtDelay = N/2;
        nCols = size(dechirpMixl,2);
        lowpassMixl = filter(digFilter,[dechirpMixl; zeros(filtDelay,nCols)]);
        lowpassMixl = lowpassMixl(filtDelay+1:end,:);
        X_3D_fullsample(:,m,l) = lowpassMixl;
    end
end
%% Downsampling
% Downsample fast-time RF samples into the ADC smaples for range FFT
dsRate = ceil(param.Fs/param.FsADC); % down sample rate
Ngds = ceil(param.FsADC*param.Tg); % number of samples per repetition period after downsampling
param.Nrfft = 2^(nextpow2(Ngds)-1); % length of range FFT
% Reduce number of samples in fast time
X_3D = complex(zeros(param.Nrfft,param.Nr,param.Nsweep));
for l = 1:param.Nsweep
    X_3D_fullsample_l = squeeze(X_3D_fullsample(:,:,l));
    X_3D_downsample_l = downsample(X_3D_fullsample_l,dsRate);
    X_3D(:,:,l) = X_3D_downsample_l(1:param.Nrfft,:); % Keep first param.Nrfft data, and discard the remaining that is typically not target samples
end
end
