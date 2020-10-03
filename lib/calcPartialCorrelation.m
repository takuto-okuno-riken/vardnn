%%
% Caluclate Partial Correlation
% returns Partial Correlation (PC) and p-values (P)
% input:
%  X     multivariate time series matrix (node x time series)

function [PC, P] = calcPartialCorrelation(X)
    [PC, P] = partialcorr(X.');
end
