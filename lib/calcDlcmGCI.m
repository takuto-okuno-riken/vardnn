%%
% Caluclate DLCM Granger Causality
% returns DLCM Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  netDLCM      trained DLCM network

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI(X, netDLCM)
    [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI_(X, [], [], [], netDLCM);
end

