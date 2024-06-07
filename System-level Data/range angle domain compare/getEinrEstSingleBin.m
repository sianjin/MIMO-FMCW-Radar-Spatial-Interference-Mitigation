% Obtain EINR estimation for a range-velocity bins from nearby
% range-velocity bins
function einr_est = getEinrEstSingleBin(n,l,hsquare_estRaw,noiseEstRaw,param,rangeGI,velocityGI,rangeTrainCells,velocityTrainCells)
% estimated h_square and noise power
hsqaure_est = zeros(param.Nt,param.numInt);
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
% covariance matrix estimation
for rangeBinIdx = 1:numExtendedRangeBins
    rangeBin = binRangeTargetExtend(rangeBinIdx);
    for velocityBinIdx = 1:numExtendedVelocityBins
        velocityBin = binVelococityTargetExtend(velocityBinIdx);
        hsqaure_est = hsqaure_est + squeeze(hsquare_estRaw(rangeBin,velocityBin,:,:));
        noise_est = noise_est + squeeze(noiseEstRaw(rangeBin,velocityBin));
    end
end
einr_est = hsqaure_est/noise_est;
end
