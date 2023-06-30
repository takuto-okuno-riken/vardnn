%%
% calculate Partial Auto Correlation of each node
% returns Partial Auto Correlation (node x PAC values)
% input:
%  X         multivariate time series matrix (node x time series)
%  maxlag    maxlag of normalized partial auto-correlation [0, maxlag] (default:15)

function [AC, lags] = calcPartialAutoCorrelation(X, maxlag)
    if nargin < 2, maxlag = 15; end
    nodeNum = size(X,1);

    AC = nan(nodeNum,maxlag+1);
%    for i=1:nodeNum
    parfor i=1:nodeNum
        x = X(i,:).';
        ulen = length(unique(single(x))); % 'half' does not support

        if ulen==1
            AC(i,:) = 0;
        else
            AC(i,:) = parcorr(double(x), NumLags=maxlag);
        end
    end
    lags = 0:maxlag;
end
