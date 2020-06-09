%%
% Caluclate Functional Connectivity
% returns Functional Connectivity (FC)
% input:
%  X     multivariate time series matrix (node x time series)

function [FC] = calcFunctionalConnectivity(X)
    FC = corr(X.',X.');
end
