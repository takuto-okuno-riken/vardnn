%%
% Estimate hidden neurons and initial weight and Create DLCM neural network
% input:
%  X          multivariate time series matrix (node x time series)
%  inSignal   multivariate time series matrix (exogenous input x time series) (optional)
%  inControl  exogenous input control matrix for each node (node x exogenous input) (optional)

function net = initDlcmNetwork(X, inSignal, inControl, bias)
    if nargin < 4, bias = 0; end
    if nargin < 3, inControl = []; end
    if nargin < 2, inSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    inputNum = size(inSignal,1);

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(nodeNum, sigLen);
    
    % set initial weight
    % For uniform distribution, bias = 0 and empty initial weight is better
    % For fMRI BOLD signal, bias = 0.5 and rough initial weight is better
    if bias > 0
        initWeight = @estimateInitWeightRoughHe;
    else
        initWeight = [];
    end
    initBias = ones(hiddenNums(1),1) * bias;

    % layer parameters
    net = createDlcmNetwork(nodeNum, inputNum, hiddenNums, inControl, initWeight, initBias);
end