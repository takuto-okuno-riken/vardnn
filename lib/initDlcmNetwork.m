%%
% Estimate hidden neurons and initial weight and Create DLCM neural network
% input:
%  X          multivariate time series matrix (node x time series)
%  inSignal   multivariate time series matrix (exogenous input x time series) (optional)
%  inControl  exogenous input control matrix for each node (node x exogenous input) (optional)

function net = initDlcmNetwork(X, inSignal, inControl)
    if nargin < 3, inControl = []; end
    if nargin < 2, inSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    inputNum = size(inSignal,1);

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(nodeNum, sigLen);
    
    % estimate init weight by granger causality
    if inputNum > 0
        X = [X; inSignal];
    end
    initWeight = estimateInitWeightByGCI(X);

    % layer parameters
    net = createDlcmNetwork(nodeNum, inputNum, hiddenNums, inControl, initWeight);
end