%%
% get VARDNN weight causality index matrix
% input:
%  net          trained VARDNN network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIv1(net, nodeControl, exControl)
    if nargin < 3
        exControl = [];
    end
    if nargin < 2
        nodeControl = [];
    end
    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;

    wcI = nan(nodeNum,nodeNum);
    for i=1:nodeNum
        % calc liner weights relation
        w1 = net.nodeNetwork{i, 1}.Layers(2, 1).Weights;

        % imparement node signals
        for j=1:nodeNum
            wcI(i,j) = var(w1(:,j),1);
        end
    end
end