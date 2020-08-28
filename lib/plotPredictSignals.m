%%
% Prot predicted signals by traind DLCM
% input:
%  X            multivariate time series matrix (node x time series)
%  inSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network
%  showEach     showSubplot (0) or show each graph (1) (default:0)

function [S, time, mae, maeerr] = plotPredictSignals(X, inSignal, nodeControl, inControl, netDLCM, showEach)
    if nargin < 6, showEach = 0; end
    [S, time, mae, maeerr] = predictDlcmNetwork(X, inSignal, nodeControl, inControl, netDLCM);
    plotTwoSignals(X, S, showEach);
end
