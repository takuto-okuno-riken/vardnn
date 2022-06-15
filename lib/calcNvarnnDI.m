%%
% Calculate nVARNN Directional Influence matrix (DI) and impaired node signals (DIsub)
% returns nVARNN DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained nVARNN network
%  nodeControl  node control matrix (1 x node) (optional)
%  exControl    exogenous input control matrix for each node (1 x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = calcNvarnnDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end

    nodeNum = net.nodeNum;
    exNum = net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end

    % set control 3D matrix (1 x node x lags)
    if isempty(nodeControl)
        nodeControl = ones(1,nodeNum,lags);
    end
    if isempty(exControl)
        exControl = ones(1,exNum,lags);
    end
    [~, ~, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc multivariate VAR DNN DI
    DI = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    nodeInputOrg = ones((nodeNum + exNum)*lags, 1);

    idx = find(control(:)==1);
    nodeInput = nodeInputOrg(idx,1);

    % predict 
    DIsub(:,1) = predict(net.network, nodeInput, 'ExecutionEnvironment', 'cpu');

    % imparement node signals
    for j=1:nodeMax
        nodeInput = nodeInputOrg;
        for p=1:lags, nodeInput(j+(nodeNum + exNum)*(p-1),:) = 0; end
        nodeInput = nodeInput(idx,1);

        % predict 
        DIsub(:,j+1) = predict(net.network, nodeInput, 'ExecutionEnvironment', 'cpu');
        DI(:,j) = abs(DIsub(:,1) - DIsub(:,j+1));
    end
end

