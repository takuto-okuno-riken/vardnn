%%
% get DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIp1(netDLCM, nodeControl, exControl)
    if nargin < 3
        exControl = [];
    end
    if nargin < 2
        nodeControl = [];
    end
    nodeNum = length(netDLCM.nodeNetwork);
    nodeInNum = size(netDLCM.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2);
    wcI = nan(nodeNum,nodeNum);
    for i=1:nodeNum
        % get input control
        control = ones(1, nodeNum);
        excontrol = ones(1, nodeInNum - nodeNum);
        if ~isempty(nodeControl)
            control = nodeControl(i,:);
        end
        if ~isempty(exControl)
            excontrol = exControl(i,:);
        end
        ctrl = [control, excontrol];

        % calc liner weights relation
        w1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        b1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Bias;
        w1 = w1 .* repmat(ctrl, size(w1,1), 1);

        v2 = w1(:,i) + b1;
        v2(v2<0) = 0; % TODO: remove?
        VarEi = var(v2(:));

        % causal inference of other nodes
        for j=1:nodeNum
            if i==j, continue; end
            v1 = w1(:,i) + w1(:,j) + b1;
            v1(v1<0) = 0; % TODO: remove?
            VarEj = var(v1(:));
            wcI(i,j) = log(VarEj / VarEi);
        end
    end
end