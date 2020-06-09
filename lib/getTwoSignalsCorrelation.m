%%
% Get correlation between two signals
% input:
%  X          multivariate time series matrix (node x time series)
%  Y          multivariate time series matrix (node x time series)

function [R] = getTwoSignalsCorrelation(X, Y)
    R = corr2(X(:), Y(:));
end
