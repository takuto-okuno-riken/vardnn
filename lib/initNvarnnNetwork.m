%%
% Estimate hidden neurons and initial weight and create Nonlinear VAR Neural Network (NVARNN)
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (1 x node) (default:[])
%  exControl       exogenous input control matrix for each node (1 x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:1)
%  activateFunc    activation function for each layer (default:@reluLayer)

function net = initNvarnnNetwork(X, exSignal, nodeControl, exControl, lags, activateFunc)
    if nargin < 6, activateFunc = @reluLayer; end % to use sigmoidLayer, we need 2020b
    if nargin < 5, lags = 1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % set control 3D matrix (1 x node x lags)
    if isempty(nodeControl)
        nodeControl = ones(1,nodeNum,lags);
    end
    if isempty(exControl)
        exControl = ones(1,exNum,lags);
    end
    [nodeControl,exControl,~] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    inputNums = sum(nodeControl,'all') + sum(exControl,'all');

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(inputNums, sigLen);
    
    % layer parameters
    net = createNvarnnNetwork(nodeNum, exNum, inputNums, hiddenNums, lags, activateFunc);
end