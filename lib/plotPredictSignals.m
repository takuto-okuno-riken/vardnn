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
    [Y, time] = predictDlcmNetwork(X, inSignal, nodeControl, inControl, netDLCM);
    S = X;
    S(:,2:end) = Y(:,1:end-1);
    A = Y(:,1:end-1) - X(:,2:end);
    mae = nanmean(abs(A(:)));
    maeerr = std(A(:),1) / sqrt(length(A(:)));
    plotTwoSignals(X, S, showEach);
end
