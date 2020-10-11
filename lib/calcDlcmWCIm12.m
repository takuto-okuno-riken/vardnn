%%
% get DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIp123(netDLCM, nodeControl, inControl)
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
        % get input control
        % calc liner weights relation
        w1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        w2 = netDLCM.nodeNetwork{i, 1}.Layers(4, 1).Weights;
        w3 = netDLCM.nodeNetwork{i, 1}.Layers(6, 1).Weights;
        b1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Bias;
        b2 = netDLCM.nodeNetwork{i, 1}.Layers(4, 1).Bias;
        b3 = netDLCM.nodeNetwork{i, 1}.Layers(6, 1).Bias;

        w = w2 * (w1 + b1);
        w(w<0) = 0; % TODO: remove?
        
        v2 = sum(w,2) + b2;
        %v2(v2<0) = 0; % TODO: remove?
        VarEi = var(v2(:));

        % imparement node signals
        for j=1:nodeNum
            if i==j, continue; end
            wt = w;
            wt(:,j) = [];
            v1 = sum(wt,2) + b2;
            %v1(v1<0) = 0; % TODO: remove?
            VarEj = var(v1(:));
            wcI(i,j) = log(VarEi / VarEj);
        end
    end
end