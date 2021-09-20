%%
% Create multivaliate VAR LSTM layers for single node
% input:
%  nodeNum        node number
%  exNum          exogenous input number
%  hiddenNums     hidden layer (next of input) neuron numbers of single unit (vector)

function layers = createMvarLstmLayers(nodeNum, exNum, hiddenNums)
    layers = [
        sequenceInputLayer(nodeNum+exNum);
        lstmLayer(hiddenNums,'OutputMode','last')
        % Add a fully connected layer
        fullyConnectedLayer(1)
        % reggression for learning
        regressionLayer();
    ];
end

