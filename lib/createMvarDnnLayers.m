%%
% Create multivaliate VAR DNN layers for single node
% input:
%  nodeNum        node number
%  exNum          exogenous input number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  lags           number of lags for autoregression
%  nNodeControl   node control matrix (1 x node) (optional)
%  nExControl     exogenous input control (1 x exogenous input) (optional)
%  activateFunc    activation function for each layer (default:@reluLayer)
%  activateFunc   activation function for each layer (optional)
%  initialWeight  weight initialize matrix of hidden1 layer (optional)
%  initialBias    bias initialize matrix of hidden1 layer (optional)

function layers = createMvarDnnLayers(nodeNum, exNum, hiddenNums, lags, nNodeControl, nExControl, activateFunc, initWeightFunc, initWeightParam, initBias, currentNode)
    if nargin < 8, initWeightFunc = []; initWeightParam = []; initBias = []; currentNode = 0; end
    if nargin < 7, activateFunc = @reluLayer; end
    if nargin < 6, nExControl = []; end
    if nargin < 5, nNodeControl = []; end

    % init first fully connected layer
    v = ver('nnet');
    nnetver = str2num(v.Version);
    if nnetver < 12.1
        firstFCLayer = fullyConnectedLayer(hiddenNums(1));
    else
        if isempty(initBias)
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializer(sz, lags, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode));
        else
            % set initial bias for each neuron
            if length(initBias) > 1
                initBias1 = ones(hiddenNums(1),1) * initBias(1);
            else
                initBias1 = ones(hiddenNums(1),1) * initBias;
            end
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializer(sz, lags, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode), ...
                'Bias', initBias1);
        end
    end
    
    %
    inLayers = [
        % input layer
        sequenceInputLayer((nodeNum+exNum)*lags);
        % Add a fully connected layer
        firstFCLayer;
        % Add an ReLU non-linearity.
        activateFunc();
        ];

    hdLayers = [];
    for i=2:length(hiddenNums)
        if nnetver < 12.1 || isempty(initBias)
            hiddenFCLayer = fullyConnectedLayer(hiddenNums(i));
        else
            if length(initBias) > 1
                initBiasH = ones(hiddenNums(i),1) * initBias(i);
            else
                initBiasH = ones(hiddenNums(i),1) * initBias;
            end
            hiddenFCLayer = fullyConnectedLayer(hiddenNums(i), 'Bias', initBiasH);
        end
        hdLayers = [
            hdLayers;
            % Add a fully connected layer
            hiddenFCLayer;
            % Add an ReLU non-linearity.
            activateFunc();
        ];
    end

    if nnetver < 12.1 || isempty(initBias)
        outputFCLayer = fullyConnectedLayer(1);
    else
        if length(initBias) > 1
            initBiasO = initBias(i+1);
        else
            initBiasO = initBias;
        end
        outputFCLayer = fullyConnectedLayer(1, 'Bias', initBiasO);
    end
    layers = [
        inLayers;
        hdLayers;
        % Add a fully connected layer
        outputFCLayer;
        % reggression for learning
        regressionLayer();
    ];
end

%%
% weight initializer
% Returns He distribution + user specified weight
function weights = weightInitializer(sz, lags, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode)
    global dlcmInitWeights;

    if ~isempty(initWeightFunc)
        weights = initWeightFunc(sz,initWeightParam);
    else
        scale = 0.1;
        filterSize = [sz(1) sz(2)];
        numIn = filterSize(1) * filterSize(2);

        varWeights = 2 / ((1 + scale^2) * numIn);
        weights = randn(sz) * sqrt(varWeights);
    end
    if ~isempty(nNodeControl)
        nodeNum = length(nNodeControl);
        filter = repmat(nNodeControl, size(weights,1), 1);
        weights(:, 1:nodeNum) = weights(:, 1:nodeNum) .* filter;
    end
    if ~isempty(nExControl)
        nodeNum = sz(2) - length(nExControl);
        filter = repmat(nExControl, size(weights,1), 1);
        weights(:, nodeNum+1:sz(2)) = weights(:, nodeNum+1:sz(2)) .* filter;
    end
    dlcmInitWeights = weights;
end
