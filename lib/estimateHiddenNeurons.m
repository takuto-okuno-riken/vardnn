%%
% Estimate neuron numbers of DLCM hidden layers
% input:
%  nodeNum    DLCM node number
%  sigLen     data time series length

function hiddenNums = estimateHiddenNeurons(nodeNum, sigLen)
    hiddenNums = zeros(2,1);
    hiddenNums(1) = ceil(32 + (sigLen-100)*0.12/(1+nodeNum*0.01));
    hiddenNums(2) = ceil(hiddenNums(1)*2/3);
%hiddenNums(3) = ceil(hiddenNums(2)*2/3);
end
