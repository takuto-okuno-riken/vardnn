%%
% get DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIv1(netDLCM, nodeControl, inControl)
    if nargin < 3
        inControl = [];
    end
    if nargin < 2
        nodeControl = [];
    end
    nodeNum = length(netDLCM.nodeNetwork);
    nodeInNum = size(netDLCM.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2);
    wcI = nan(nodeNum,nodeNum);
    for i=1:nodeNum
        % calc liner weights relation
        w1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;

        % imparement node signals
        for j=1:nodeNum
            wcI(i,j) = var(w1(:,j));
        end
    end
end