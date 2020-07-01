%%
% Estimate hidden neurons and initial weight and Create DLCM neural network
% input:
%  X               multivariate time series matrix (node x time series)
%  inSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  inControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  initWeightFunc  initializing weight function (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = initDlcmNetwork(X, inSignal, inControl, initWeightFunc, initWeightParam, initBias)
    if nargin < 6, initBias = 0; end
    if nargin < 5, initWeightParam = []; end
    if nargin < 4, initWeightFunc = []; end
    if nargin < 3, inControl = []; end
    if nargin < 2, inSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    inputNum = size(inSignal,1);

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(nodeNum, sigLen);
    
    % set initial bias for each neuron
    biasMat = ones(hiddenNums(1),1) * initBias;

    % layer parameters
    net = createDlcmNetwork(nodeNum, inputNum, hiddenNums, inControl, initWeightFunc, initWeightParam, biasMat);
end