%%
% Create DLCM neural network
% input:
%  nodeNum         DLCM node number
%  exNum           DLCM exogenous input number
%  hiddenNums      hidden layer (next of input) neuron numbers of single unit (vector)
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix (node x exogenous input) (default:[])
%  initWeightFunc  initializing weight function for hidden1 layer (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value for hidden1 layer (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = createDlcmNetwork(nodeNum, exNum, hiddenNums, nodeControl, exControl, initWeightFunc, initWeightParam, initBias)
    if nargin < 8, initBias = zeros(hiddenNums(1),1); end
    if nargin < 7, initWeightParam = []; end
    if nargin < 6, initWeightFunc = []; end
    if nargin < 5, exControl = []; end
    if nargin < 4, nodeControl = []; end

    nodeLayers = cell(nodeNum,1);
    for i=1:nodeNum
        nNodeControl = [];
        if ~isempty(nodeControl), nNodeControl = nodeControl(i,:); end
        nExControl = [];
        if ~isempty(exControl), nExControl = exControl(i,:); end
        nodeLayers{i} = createDlcmLayers(nodeNum, exNum, hiddenNums, nNodeControl, nExControl, initWeightFunc, initWeightParam, initBias, i);
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,1);
    net.trainInfo = cell(nodeNum,1);
    net.initWeights = cell(nodeNum,1);
    net.trainTime = 0;
    net.recoverTrainTime = 0;
end