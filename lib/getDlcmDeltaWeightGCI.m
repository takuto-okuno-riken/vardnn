%%
% get DLCM delta weight granger causality index matrix
% input:
%  netDLCM     trained DLCM network

function gcI = getDlcmDeltaWeightGCI(netDLCM)
    nodeNum = netDLCM.nodeNum;
    nodeInNum = nodeNum + net.exNum;

    gcI = nan(nodeNum,nodeInNum);
    for i=1:nodeNum
        weights = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        initWeights = netDLCM.initWeights{i};
        weight = weights - initWeights;
        VarEi = var(weight(:));

        % imparement node signals
        for j=1:nodeInNum
            if i==j, continue; end
            nweight = weight;
            nweight(:,j) = [];
            VarEj = var(nweight(:));
            gcI(i,j) = log(VarEi / VarEj);
        end
    end
end