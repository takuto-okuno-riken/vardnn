%%
% Estimate neuron numbers of VARDNN hidden layers
% input:
%  nodeNum         VARDNN node number
%  sigLen          data time series length
%  maxNeuronNum    maximum neuron number of a hidden layer (default:NaN)

function hiddenNums = estimateHiddenNeurons(nodeNum, sigLen, maxNeuronNum)
    if nargin < 3, maxNeuronNum = nan; end
    hiddenNums = zeros(2,1);
    hiddenNums(1) = ceil(32 + (sigLen-100)*0.12/(1+nodeNum*0.01));
    if maxNeuronNum > 0 && hiddenNums(1) > maxNeuronNum, hiddenNums(1) = maxNeuronNum; end
    hiddenNums(2) = ceil(hiddenNums(1)*2/3);
%hiddenNums(3) = ceil(hiddenNums(2)*2/3);
end
