%%
% Create DLCM neural network layers for single node
% input:
%  nodeNum        DLCM node number
%  inputNum       DLCM exogenous input number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  nodeInControl  exogenous input control (1 x exogenous input) (optional)
%  initialWeight  weight initialize matrix of hidden1 layer (optional)
%  initialBias    bias initialize matrix of hidden1 layer (optional)

function layers = createDlcmLayers(nodeNum, inputNum, hiddenNums, nodeInControl, initialWeight, initBias, currentNode)
    if nargin < 6, initialWeight = []; currentNode = 0; end
    if nargin < 5, nodeInControl = []; end

    % init first fully connected layer
    if isempty(initialWeight) && isempty(nodeInControl) && isempty(initBias)
        firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
            'WeightsInitializer', @(sz) weightedHe(sz, nodeInControl, initialWeight, currentNode));
    else
        firstFCLayer = fullyConnectedLayer(hiddenNums(1), ...
            'WeightsInitializer', @(sz) weightedHe(sz, nodeInControl, initialWeight, currentNode), ...
            'Bias', initBias);
    end
    
    %
    inLayers = [
        % input layer
        sequenceInputLayer(nodeNum+inputNum);
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
    ]
end

%%
% weight initializer
% Returns He distribution + user specified weight
function weights = weightedHe(sz, nodeInControl, initialWeight, currentNode)
    global dlcmInitWeights;

    scale = 0.1;

    filterSize = [sz(1) sz(2)];
    numIn = filterSize(1) * filterSize(2);

    varWeights = 2 / ((1 + scale^2) * numIn);
    weights = randn(sz) * sqrt(varWeights);
    if ~isempty(nodeInControl)
        nodeNum = sz(2) - length(nodeInControl);
        filter = repmat(nodeInControl, size(weights,1), 1);
        weights(:, nodeNum+1:sz(2)) = weights(:, nodeNum+1:sz(2)) .* filter;
    end
    if ~isempty(initialWeight)
        A = initialWeight(currentNode,:);
        if ~isempty(nodeInControl)
            nodeNum = sz(2) - length(nodeInControl);
            A(:, nodeNum+1:sz(2)) = A(:, nodeNum+1:sz(2)) .* nodeInControl;
        end
        B = repmat(A,sz(1),1);
        weights = weights + B;
    end
    dlcmInitWeights = weights;
end
