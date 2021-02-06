%%
% Return trained DLCM network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           DLCM network structure
%  options       training options

function trainedNet = trainDlcmNetwork(X, exSignal, nodeControl, exControl, net, options)
    trainedNet = trainMvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options);
end