%%
% get VARDNN weight causality index matrix
% input:
%  net          trained VARDNN network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCIm12(net, nodeControl, exControl)
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
        w1 = net.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        w2 = net.nodeNetwork{i, 1}.Layers(4, 1).Weights;
        w3 = net.nodeNetwork{i, 1}.Layers(6, 1).Weights;
        b1 = net.nodeNetwork{i, 1}.Layers(2, 1).Bias;
        b2 = net.nodeNetwork{i, 1}.Layers(4, 1).Bias;
        b3 = net.nodeNetwork{i, 1}.Layers(6, 1).Bias;
        w1 = w1 .* repmat(ctrl, size(w1,1), 1);

        x = sum(w1,2) + b1;
        x(x<0) = 0;

        y = w2 * x + b2;
        y(y<0) = 0;
        VarEi = var(y(:));
        
        % imparement node signals
        for j=1:nodeNum
            if i==j, continue; end
            w = w1;
            w(:,j) = 0;
            xj = sum(w,2) + b1;
            xj(xj<0) = 0;
            
            yj = w2 * xj + b2;
            yj(yj<0) = 0;
            VarEj = var(yj(:));
            wcI(i,j) = log(VarEi / VarEj);
        end
    end
end