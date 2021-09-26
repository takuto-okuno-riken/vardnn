%%
% Estimate hidden neurons and initial weight and create pairwise VAR DNN
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  initWeightFunc  initializing weight function (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = initPvarDnnNetwork(X, exSignal, nodeControl, exControl, lags, initWeightFunc, initWeightParam, initBias)
    if nargin < 8, initBias = 0; end
    if nargin < 7, initWeightParam = []; end
    if nargin < 6, initWeightFunc = []; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,~] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    inputNums = ceil(sum(nodeControl,'all') / (nodeNum * nodeNum)) * 2;

    % estimate neuron number of hidden layers
    hiddenNums = zeros(2,1);
    hiddenNums(1) = ceil(16 + (sigLen-100)*0.12/(1+inputNums*0.01));
    hiddenNums(2) = ceil(hiddenNums(1)*2/3);
    
    % set initial bias for each neuron
    biasMat = ones(hiddenNums(1),1) * initBias;

    % layer parameters
    net = createPvarDnnNetwork(nodeNum, exNum, hiddenNums, lags, nodeControl, exControl, initWeightFunc, initWeightParam, biasMat);
end