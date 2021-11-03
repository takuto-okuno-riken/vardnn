%%
% Create Nonlinear VAR Neural Network (NVARNN) network
% input:
%  nodeNum         node number
%  exNum           exogenous input number
%  inputNums       total input number (through control)
%  hiddenNums      hidden layer (next of input) neuron numbers of single unit (vector)
%  lags            number of lags for autoregression (default:1)
%  activateFunc    activation function for each layer (default:@reluLayer)

function net = createNvarnnNetwork(nodeNum, exNum, inputNums, hiddenNums, lags, activateFunc)
    if nargin < 6, activateFunc = @reluLayer; end % to use sigmoidLayer, we need 2020b
    if nargin < 5, lags = 1; end

    % create layers
    layers = createNvarnnLayers(nodeNum, exNum, inputNums, hiddenNums, activateFunc);

    % set structure
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.layers = layers;
    net.network = [];
    net.trainInfo = [];
    net.initWeights = [];
    net.trainTime = 0;
end