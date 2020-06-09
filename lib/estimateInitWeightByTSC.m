%%
% Estimate initial weight for DLCM by Characterized Sliding Correlation
% input:
%  X     multivariate time series matrix (node x time series)
%  p     number of lags for autoregression (default:3)

function initWeight = estimateInitWeightByTSC(X, p)
    if nargin < 2
        p = 3;
    end
    initWeight = calcTimeShiftedCorrelation(X, p); % return is [-1 1]
    initWeight = initWeight * 0.2; % set appropriate value for initial weight
end