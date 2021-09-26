%%
% Create pairwise VAR DNN layers for single node
% input:
%  nodeNum        pairwised VAR DNN node number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  lags           number of lags for autoregression (default:3)
%  initialWeight  weight initialize matrix of hidden1 layer (optional)
%  initialBias    bias initialize matrix of hidden1 layer (optional)

function layers = createPvarDnnLayers(inputNums, hiddenNums, initWeightFunc, initWeightParam, initBias)
    if nargin < 5, initWeightFunc = []; initWeightParam = []; initBias = []; end

    % init first fully connected layer
    v = ver('nnet');
    nnetver = str2num(v.Version);
    if nnetver < 12.1
        firstFCLayer = fullyConnectedLayer(hiddenNums(1));
    else
        if isempty(initWeightFunc) && isempty(initBias)
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializerPvarDnn(sz, initWeightFunc, initWeightParam));
        else
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializerPvarDnn(sz, initWeightFunc, initWeightParam), ...
                'Bias', initBias);
        end
    end
    
    %
    inLayers = [
        % input layer
        sequenceInputLayer(inputNums);
        % Add a fully connected layer
        firstFCLayer;
        % Add an ReLU non-linearity.
        reluLayer();
        ];

    hdLayers = [];
    for i=2:length(hiddenNums)
        hdLayers = [
            hdLayers;
            % Add a fully connected layer
            fullyConnectedLayer(hiddenNums(i));
            % Add an ReLU non-linearity.
            reluLayer();
        ];
    end

    layers = [
        inLayers;
        hdLayers;
        % Add a fully connected layer
        fullyConnectedLayer(1);

        % Add an ReLU non-linearity.
        % reluLayer();

        % reggression for learning
        regressionLayer();
    ];
end

%%
% weight initializer
% Returns He distribution + user specified weight
function weights = weightInitializerPvarDnn(sz, initWeightFunc, initWeightParam)
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
