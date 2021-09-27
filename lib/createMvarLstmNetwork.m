%%
% Create multivaliate VAR LSTM network
% input:
%  nodeNum         node number
%  exNum           exogenous input number
%  hiddenNums      hidden layer (next of input) neuron numbers of single unit (vector)
%  lags            number of lags for autoregression (default:1)

function net = createMvarLstmNetwork(nodeNum, exNum, hiddenNums, lags)
    if nargin < 4, lags = 1; end

    nodeLayers = cell(nodeNum,1);
    for i=1:nodeNum
        nodeLayers{i} = createMvarLstmLayers(nodeNum, exNum, hiddenNums);
    end
    net.version = 1.1;
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,1);
    net.trainInfo = cell(nodeNum,1);
    net.trainTime = 0;
end