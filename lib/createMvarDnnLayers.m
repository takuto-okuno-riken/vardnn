%%
% Create multivaliate VAR DNN layers for single node
% input:
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  nNodeControl   node control matrix (1 x node)
%  nExControl     exogenous input control (1 x exogenous input)
%  activateFunc   activation function for each layer (default:@reluLayer)
%  initialWeight  weight initialize matrix of hidden1 layer (optional)
%  initialBias    bias initialize matrix of hidden1 layer (optional)

function layers = createMvarDnnLayers(hiddenNums, nNodeControl, nExControl, activateFunc, initWeightFunc, initWeightParam, initBias, currentNode)
    if nargin < 5, initWeightFunc = []; initWeightParam = []; initBias = []; currentNode = 0; end
    if nargin < 4, activateFunc = @reluLayer; end

    % init first fully connected layer
    v = ver('nnet');
    nnetver = str2num(v.Version);
    if nnetver < 12.1
        firstFCLayer = fullyConnectedLayer(hiddenNums(1));
    else
        if isempty(initBias)
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializer(sz, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode));
        else
            % set initial bias for each neuron
            if length(initBias) > 1
                initBias1 = ones(hiddenNums(1),1) * initBias(1);
            else
                initBias1 = ones(hiddenNums(1),1) * initBias;
            end
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializer(sz, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode), ...
                'Bias', initBias1);
        end
    end
    
    %
    inputNums = sum(nNodeControl,'all') + sum(nExControl,'all');
    inLayers = [
        % input layer
        sequenceInputLayer(inputNums);
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
function weights = weightInitializer(sz, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode)
    global dnnInitWeights;

    if ~isempty(initWeightFunc)
        weights = initWeightFunc(sz,initWeightParam);
    else
        scale = 0.1;
        filterSize = [sz(1) sz(2)];
        numIn = filterSize(1) * filterSize(2);

        varWeights = 2 / ((1 + scale^2) * numIn);
        weights = randn(sz) * sqrt(varWeights);
    end
    dnnInitWeights = weights;
end
