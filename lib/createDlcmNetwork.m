%%
% Create DLCM neural network
% input:
%  nodeNum        DLCM node number
%  inputNum       DLCM exogenous input number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  inControl      exogenous input control matrix (node x exogenous input) (optional)
%  initialWeight  weight initialize matrix of hidden1 layer (optional)

function net = createDlcmNetwork(nodeNum, inputNum, hiddenNums, inControl, initialWeight)
    if nargin < 6, initialWeight = []; end
    if nargin < 5, inControl = []; end

    nodeLayers = cell(nodeNum,1);
    for i=1:nodeNum
        nodeInControl = [];
        if ~isempty(inControl), nodeInControl = inControl(i,:); end
        nodeLayers{i} = createDlcmLayers(nodeNum, inputNum, hiddenNums, nodeInControl, initialWeight, i);
    end
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,1);
    net.trainInfo = cell(nodeNum,1);
    net.initWeights = cell(nodeNum,1);
    net.trainTime = 0;
    net.recoverTrainTime = 0;
end