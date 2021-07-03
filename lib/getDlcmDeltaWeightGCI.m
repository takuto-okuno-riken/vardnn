%%
% get VARDNN delta weight granger causality index matrix
% input:
%  net         trained VARDNN network

function gcI = getDlcmDeltaWeightGCI(net)
    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;

    gcI = nan(nodeNum,nodeInNum);
    for i=1:nodeNum
        weights = net.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        initWeights = net.initWeights{i};
        weight = weights - initWeights;
        VarEi = var(weight(:),1);

        % imparement node signals
        for j=1:nodeInNum
            if i==j, continue; end
            nweight = weight;
            nweight(:,j) = [];
            VarEj = var(nweight(:),1);
            gcI(i,j) = log(VarEi / VarEj);
        end
    end
end