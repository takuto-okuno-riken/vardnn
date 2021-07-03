%%
% get VARDNN DI matrix (vdDI) and its standard error matrix
% input:
%  net     trained VARDNN network

function [vdDI] = getDlcmECcorrDeltaWeight(net)
    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;

    vdDI = zeros(nodeNum,nodeInNum);
    for i=1:nodeNum
        weights = net.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        initWeights = net.initWeights{i};
        for j=1:nodeInNum
            vdDI(i,j) = getCosSimilarity(weights(:,j),initWeights(:,j)); % corr2
        end
    end
end