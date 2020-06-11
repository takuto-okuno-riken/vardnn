%%
% get DLCM weight granger causality index matrix
% input:
%  netDLCM     trained DLCM network

function gcI = getDlcmWeightGCI(netDLCM)
    nodeNum = length(netDLCM.nodeNetwork);
    nodeInNum = size(netDLCM.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2);
    gcI = nan(nodeNum,nodeInNum);
    for i=1:nodeNum
        weight = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
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