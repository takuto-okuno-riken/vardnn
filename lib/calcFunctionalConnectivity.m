%%
% Caluclate Functional Connectivity
% returns Functional Connectivity (FC) and p-values (P)
% input:
%  X     multivariate time series matrix (node x time series)

function [FC, P] = calcFunctionalConnectivity(X)
    [FC, P] = corr(X.', X.');
end
