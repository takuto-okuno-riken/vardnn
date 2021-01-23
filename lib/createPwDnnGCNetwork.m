%%
% Create Pairwised DNN-GC's neural network
% input:
%  nodeNum         node number
%  exNum           exogenous input number
%  hiddenNums      hidden layer (next of input) neuron numbers of single unit (vector)
%  lags            number of lags for autoregression (default:3)
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix (node x exogenous input) (default:[])
%  initWeightFunc  initializing weight function for hidden1 layer (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value for hidden1 layer (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = createPwDnnGCNetwork(nodeNum, exNum, hiddenNums, lags, nodeControl, exControl, initWeightFunc, initWeightParam, initBias)
    if nargin < 9, initBias = zeros(hiddenNums(1),1); end
    if nargin < 8, initWeightParam = []; end
    if nargin < 7, initWeightFunc = []; end
    if nargin < 6, exControl = []; end
    if nargin < 5, nodeControl = []; end
    if nargin < 4, lags = 3; end
    
    nodeMax = nodeNum + exNum;
    nodeLayers = cell(nodeNum,nodeMax);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            nodeLayers{i,j} = createPwDnnGCLayers(2, hiddenNums, lags, initWeightFunc, initWeightParam, initBias);
        end
    end
    net.lags = lags;
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,nodeMax);
    net.trainInfo = cell(nodeNum,nodeMax);
    net.initWeights = cell(nodeNum,nodeMax);
    net.trainTime = 0;
    net.recoverTrainTime = 0;
end