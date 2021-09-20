%%
% Estimate hidden neurons and initial weight and create multivariate VARLSTM
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:1)

function net = initMvarLstmNetwork(X, exSignal, nodeControl, exControl, lags)
    if nargin < 5, lags = 1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons((nodeNum+exNum)*lags, sigLen);
    
    % layer parameters
    net = createMvarLstmNetwork(nodeNum, exNum, hiddenNums(1), lags);
end