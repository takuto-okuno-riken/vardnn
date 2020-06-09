%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM     trained DLCM network

function [dlEC, dlECerr] = getDlcmECmeanDeltaWeight(netDLCM)
    nodeNum = length(netDLCM.nodeNetwork);
    nodeInNum = size(netDLCM.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2);
    dlEC = zeros(nodeNum,nodeInNum);
    dlECerr = zeros(nodeNum,nodeInNum);
    for i=1:nodeNum
        weights = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        initWeights = netDLCM.initWeights{i};
        mweights = mean(weights-initWeights,1);
        eweights = std(weights,1) / sqrt(size(weights,1));
        dlEC(i,:) = mweights;
        dlECerr(i,:) = eweights;
    end
end