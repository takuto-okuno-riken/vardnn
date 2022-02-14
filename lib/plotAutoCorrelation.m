%%
% Plot Auto Correlation of each node
% returns Auto Correlation (node x AC values)
% input:
%  X         multivariate time series matrix (node x time series)
%  maxlag    maxlag of normalized auto-correlation [0, maxlag] (default:15)

function [AC, lags] = plotAutoCorrelation(X, maxlag)
    if nargin < 2, maxlag = 15; end

    [AC, lags] = calcAutoCorrelation(X, maxlag);

    % plot Auto Correlation
    plot(lags, AC.');
    title('Auto Correlation')
    xlabel('lags')
end
