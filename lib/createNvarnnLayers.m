%%
% Create Nonlinear VAR Neural Network (NVARNN) layers for single node
% input:
%  nodeNum        node number
%  exNum          exogenous input number
%  inputNums      total input number (through control)
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)
%  activateFunc   activation function for each layer

function layers = createNvarnnLayers(nodeNum, exNum, inputNums, hiddenNums, activateFunc)
    layers = [
        % input layer
        sequenceInputLayer(inputNums);
        % Add a fully connected layer
        fullyConnectedLayer(hiddenNums(1));
        % Add an sigmoid non-linearity.
        activateFunc();

        % Add a fully connected layer
        fullyConnectedLayer(nodeNum)
        % reggression for learning
        regressionLayer();
    ];
end

