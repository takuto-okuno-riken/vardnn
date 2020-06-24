%%
% Create DLCM neural network
% input:
%  nodeNum        DLCM node number
%  inputNum       DLCM exogenous input number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  inControl      exogenous input control matrix (node x exogenous input) (default:[])
%  initWeightFunc initializing weight function for hidden1 layer (default:[])
%  initBias       initializing bias value for hidden1 layer (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = createDlcmNetwork(nodeNum, inputNum, hiddenNums, inControl, initialWeightFunc, initBias)
    if nargin < 6, initBias = zeros(hiddenNums(1),1); end
    if nargin < 5, initialWeightFunc = []; end
    if nargin < 4, inControl = []; end

    nodeLayers = cell(nodeNum,1);
    for i=1:nodeNum
        nodeInControl = [];
        if ~isempty(inControl), nodeInControl = inControl(i,:); end
        nodeLayers{i} = createDlcmLayers(nodeNum, inputNum, hiddenNums, nodeInControl, initialWeightFunc, initBias, i);
    end
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,1);
    net.trainInfo = cell(nodeNum,1);
    net.initWeights = cell(nodeNum,1);
    net.trainTime = 0;
    net.recoverTrainTime = 0;
end