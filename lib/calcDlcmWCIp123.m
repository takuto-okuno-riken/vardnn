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
%{        
        control = ones(1, nodeNum);
        excontrol = ones(1, nodeInNum - nodeNum);
        if ~isempty(nodeControl)
            control = nodeControl(i,:);
        end
        if ~isempty(inControl)
            excontrol = inControl(i,:);
        end
        ctrl = [control, excontrol];
%}
        % calc liner weights relation
        w1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        w2 = netDLCM.nodeNetwork{i, 1}.Layers(4, 1).Weights;
        w3 = netDLCM.nodeNetwork{i, 1}.Layers(6, 1).Weights;
        b1 = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Bias;
        b2 = netDLCM.nodeNetwork{i, 1}.Layers(4, 1).Bias;
        b3 = netDLCM.nodeNetwork{i, 1}.Layers(6, 1).Bias;

        w23 = w2.' * w3.';
        b23 = b2.' * w3.' + b3;
        
%        w123 = w1.' * w23;
%        b123 = b1.' * w23 + b23;
%        w = w123.';
%        b = b123;

        w23r = repmat(w23, 1, size(w1,2));
        w = w1 .* w23r;
        b = b1 .* w23  + b23;
        
        v2 = w(:,i) + b;
        VarEi = var(v2(:));

        % imparement node signals
        for j=1:nodeNum
            if i==j, continue; end
            v1 = w(:,i) + w(:,j) + b;
            VarEj = var(v1(:));
            wcI(i,j) = log(VarEj / VarEi);
        end
    end
end