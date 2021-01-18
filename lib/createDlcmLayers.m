%%
% Create DLCM neural network layers for single node
% input:
%  nodeNum        DLCM node number
%  exNum          DLCM exogenous input number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  nNodeControl   node control matrix (1 x node) (optional)
%  nExControl     exogenous input control (1 x exogenous input) (optional)
%  initialWeight  weight initialize matrix of hidden1 layer (optional)
%  initialBias    bias initialize matrix of hidden1 layer (optional)

function layers = createDlcmLayers(nodeNum, exNum, hiddenNums, nNodeControl, nExControl, initWeightFunc, initWeightParam, initBias, currentNode)
    if nargin < 7, initWeightFunc = []; initWeightParam = []; initBias = []; currentNode = 0; end
    if nargin < 6, nExControl = []; end
    if nargin < 5, nNodeControl = []; end

    % init first fully connected layer
    v = ver('nnet');
    nnetver = str2num(v.Version);
    if nnetver < 12.1
        firstFCLayer = fullyConnectedLayer(hiddenNums(1));
    else
        if isempty(initWeightFunc) && isempty(nExControl) && isempty(initBias)
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializer(sz, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode));
        else
            firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
                'WeightsInitializer', @(sz) weightInitializer(sz, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode), ...
                'Bias', initBias);
        end
    end
    
    %
    inLayers = [
        % input layer
        sequenceInputLayer(nodeNum+exNum);
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
function weights = weightInitializer(sz, nNodeControl, nExControl, initWeightFunc, initWeightParam, currentNode)
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
