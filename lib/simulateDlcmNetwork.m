%%
% Simulate node signals by traind DLCM and exogenous input
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network

function [S, time] = simulateDlcmNetwork(X, exSignal, nodeControl, exControl, netDLCM)
    [S, time] = simulateMvarDnnNetwork(X, exSignal, nodeControl, exControl, netDLCM);
end
