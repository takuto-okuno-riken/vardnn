%%
% fitting DNN regression
% input:
%  X               regression explanatory variable (data point x expNum)
%  Y               regression objective variable (data point x 1)
%  options         training options
%  activateFunc    activation function for each layer (default:@reluLayer)

function mdl = fitDnnRegression(X, Y, options, activateFunc)
    if nargin < 4, activateFunc = @reluLayer; end
    predNum = size(X,2);

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(predNum, 0);
    layers = createDnnRegressionLayers(predNum, hiddenNums, activateFunc);
    
    % train deep neural network
    [network, trainInfo] = trainNetwork(X.', Y.', layers, options);
    
    mdl.layers = layers;
    mdl.network = network;
    mdl.trainInfo = trainInfo;
end