%%
% Create DNN Regression layers for single node
% input:
%  inputNums      input neuron numbers
%  hiddenNums     hidden layer (next of input) neuron numbers (vector)
%  activateFunc   activation function for each layer (default:@reluLayer)

function layers = createDnnRegressionLayers(inputNums, hiddenNums, activateFunc)
    % init first fully connected layer
    firstFCLayer = fullyConnectedLayer(hiddenNums(1));

    %
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
        hiddenFCLayer = fullyConnectedLayer(hiddenNums(i));
        hdLayers = [
            hdLayers;
            % Add a fully connected layer
            hiddenFCLayer;
            % Add an ReLU non-linearity.
            activateFunc();
        ];
    end

    outputFCLayer = fullyConnectedLayer(1);
    layers = [
        inLayers;
        hdLayers;
        % Add a fully connected layer
        outputFCLayer;
        % reggression for learning
        regressionLayer();
    ];
end
