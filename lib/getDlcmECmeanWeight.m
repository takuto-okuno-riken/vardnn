%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM     trained DLCM network

function [dlEC, dlECerr] = getDlcmECmeanWeight(netDLCM)
    nodeNum = netDLCM.nodeNum;
    nodeInNum = nodeNum + net.exNum;

    dlEC = zeros(nodeNum,nodeInNum);
    dlECerr = zeros(nodeNum,nodeInNum);
    for i=1:nodeNum
        weight = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        mweight = mean(weight,1);
        eweight = std(weight,1) / sqrt(size(weight,1));
        dlEC(i,:) = mweight;
        dlECerr(i,:) = eweight;
    end
end