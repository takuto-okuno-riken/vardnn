%%
% Create multivaliate VAR DNN network
% input:
%  nodeNum         node number
%  exNum           exogenous input number
%  hiddenNums      hidden layer (next of input) neuron numbers of single unit (vector)
%  lags            number of lags for autoregression (default:1)
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix (node x exogenous input) (default:[])
%  activateFunc    activation function for each layer (default:@reluLayer)
%  initWeightFunc  initializing weight function for hidden1 layer (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value for hidden1 layer (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = createMvarDnnNetwork(nodeNum, exNum, hiddenNums, lags, nodeControl, exControl, activateFunc, initWeightFunc, initWeightParam, initBias)
    if nargin < 10, initBias = zeros(hiddenNums(1),1); end
    if nargin < 9, initWeightParam = []; end
    if nargin < 8, initWeightFunc = []; end
    if nargin < 7, activateFunc = @reluLayer; end
    if nargin < 6, exControl = []; end
    if nargin < 5, nodeControl = []; end
    if nargin < 4, lags = 1; end

    nodeLayers = cell(nodeNum,1);
    for i=1:nodeNum
        nNodeControl = [];
        if ~isempty(nodeControl), nNodeControl = nodeControl(i,:); end
        nExControl = [];
        if ~isempty(exControl), nExControl = exControl(i,:); end
        nodeLayers{i} = createMvarDnnLayers(nodeNum, exNum, hiddenNums, lags, nNodeControl, nExControl, activateFunc, initWeightFunc, initWeightParam, initBias, i);
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,1);
    net.trainInfo = cell(nodeNum,1);
    net.initWeights = cell(nodeNum,1);
    net.trainTime = 0;
    net.recoverTrainTime = 0;
end