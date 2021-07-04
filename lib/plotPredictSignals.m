%%
% Prot predicted signals by traind VARDNN
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained mVAR DNN network
%  showEach     showSubplot (0) or show each graph (1) (default:0)

function [S, time, mae, maeerr] = plotPredictSignals(X, exSignal, nodeControl, exControl, net, showEach)
    if nargin < 6, showEach = 0; end
    [Y, time] = predictMvarDnnNetwork(X, exSignal, nodeControl, exControl, net);
    S = X;
    S(:,2:end) = Y(:,1:end-1);
    A = Y(:,1:end-1) - X(:,2:end);
    mae = nanmean(abs(A(:)));
    maeerr = std(A(:),1) / sqrt(length(A(:)));
    plotTwoSignals(X, S, showEach);
end
