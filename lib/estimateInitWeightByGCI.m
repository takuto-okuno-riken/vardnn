%%
% Estimate initial weight for VARDNN by Granger Causality Index (no more use this)
% input:
%  X      multivariate time series matrix (node x time series)
%  p      number of lags for autoregression (default:3)
%  range  threshold of index range (default:10)

function initWeight = estimateInitWeightByGCI(X, p, range)
    if nargin < 3
        range = 10;
    end
    if nargin < 2
        p = 3;
    end
    initWeight = calcPairwiseGCI(X, [], [], [], p);
    idx = find(~isfinite(initWeight)); % remove inf or nan
    initWeight(idx) = 0;
    idx = find(initWeight>=range); % set maxmum
    initWeight(idx) = range;
    idx = find(initWeight<=-range); % set minimum
    initWeight(idx) = -range;
    initWeight = initWeight / (range*5); % set appropriate value for initial weight
end