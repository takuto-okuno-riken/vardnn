%%
% Caluclate DLCM Granger Causality
% returns DLCM Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI(X, exSignal, nodeControl, exControl, netDLCM, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end
    [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMvarDnnGCI(X, exSignal, nodeControl, exControl, netDLCM, alpha, isFullNode);
end

