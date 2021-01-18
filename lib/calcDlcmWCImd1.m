%%
% get DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIm1(netDLCM, nodeControl, exControl)
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

        % remove useless weights
        w = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        bias = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Bias;
        for j=nodeInNum:-1:1
            if ctrl(j) < 1
                w(:,j) = [];
            end
        end
        weight = w; %[w, bias];
        VarEi = var(weight(:));

        % imparement node signals
        for j=1:nodeNum
            if i==j, continue; end
            nweight = weight;
            nweight(:,j) = [];
            VarEj = var(nweight(:));
            wcI(i,j) = log(VarEi / VarEj);
        end
    end
end