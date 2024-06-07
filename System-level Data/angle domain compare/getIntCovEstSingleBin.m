function [R_est_normalize,noise_est] = getIntCovEstSingleBin(Y_rD_3D,n,l,noiseEstRaw,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells,param)
% estimated noise power
noise_est = 0;
% range bins for estimation
binRangeTargetExtend = [max(n-rangeGI-rangeTrainCells,1):n-1-rangeGI, ...
    n+1+rangeGI:min(n+rangeGI+rangeTrainCells,size(noiseEstRaw,1))]; % range bins near the target
% velocity bins for estimation
binVelococityTargetExtend = [max(l-velocityGI-velocityTrainCells,1):l-1-velocityGI, ...
    l+1+velocityGI:min(l+velocityGI+velocityTrainCells,size(noiseEstRaw,2))]; % velocity bins near the target
% number of range bins for estimation
numExtendedRangeBins = length(binRangeTargetExtend);
% number of velocity bins for estimation
numExtendedVelocityBins = length(binVelococityTargetExtend);
% Noise estimation
for rangeBinIdx = 1:numExtendedRangeBins
    rangeBin = binRangeTargetExtend(rangeBinIdx);
    for velocityBinIdx = 1:numExtendedVelocityBins
        velocityBin = binVelococityTargetExtend(velocityBinIdx);
        noise_est = noise_est + squeeze(noiseEstRaw(rangeBin,velocityBin));
    end
end
noise_est = noise_est/numExtendedRangeBins/numExtendedVelocityBins;
% Covariance matrix estimation for range bin n and velocity bin l
R_est = zeros(param.Nt*param.Nr);
for rangeBinIdx = 1:numExtendedRangeBins
    rangeBin = binRangeTargetExtend(rangeBinIdx);
    for velocityBinIdx = 1:numExtendedVelocityBins
        velocityBin = binVelococityTargetExtend(velocityBinIdx);
        y_n_l = squeeze(Y_rD_3D(rangeBin,:,velocityBin))';
        R_est = R_est + y_n_l*y_n_l';
    end
end
R_est = R_est/numExtendedRangeBins/numExtendedVelocityBins;
R_est_normalize = R_est/noise_est;
end
