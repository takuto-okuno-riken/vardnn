%%
% get DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIm1(netDLCM, nodeControl, inControl)
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
        b1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Bias;

        v2 = sum(w1,2) + b1;
        v2(v2<0) = 0; % TODO: remove?
        VarEi = var(v2(:));

        % causal inference of other nodes
        for j=1:nodeNum
            if i==j, continue; end
            w = w1;
            w(:,j) = [];
            v1 = sum(w,2) + b1;
            v1(v1<0) = 0; % TODO: remove?
            VarEj = var(v1(:));
            wcI(i,j) = log(VarEi / VarEj);
        end
    end
end