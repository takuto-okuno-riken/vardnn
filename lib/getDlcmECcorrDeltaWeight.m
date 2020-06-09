%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM     trained DLCM network

function [dlEC] = getDlcmECcorrDeltaWeight(netDLCM)
    nodeNum = length(netDLCM.nodeNetwork);
    nodeInNum = size(netDLCM.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2);
    dlEC = zeros(nodeNum,nodeInNum);
    for i=1:nodeNum
        weights = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        initWeights = netDLCM.initWeights{i};
        for j=1:nodeInNum
            dlEC(i,j) = getCosSimilarity(weights(:,j),initWeights(:,j)); % corr2
        end
    end
end