%%
% calculate Auto Correlation of each node
% returns Auto Correlation (node x AC values)
% input:
%  X         multivariate time series matrix (node x time series)
%  maxlag    maxlag of normalized auto-correlation [0, maxlag] (default:15)

function [AC, lags] = calcAutoCorrelation(X, maxlag)
    if nargin < 2, maxlag = 15; end
    nodeNum = size(X,1);

    AC = nan(nodeNum,maxlag+1);
    for i=1:nodeNum
        x = X(i,:).';
        [c, lags] = xcov(x,x,maxlag,'normalized');
        AC(i,:) = c(maxlag+1:end);
    end
    lags = lags(maxlag+1:end);
end
