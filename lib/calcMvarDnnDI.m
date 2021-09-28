%%
% Caluclate multivariate VAR DNN Directional Influence matrix (DI) and impaired node signals (DIsub)
% returns mVAR DNN DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained multivariate VAR DNN network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = calcMvarDnnDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end

    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    if isfield(net, 'exNum'), exNum = net.exNum; else exNum = size(net.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2) - nodeNum; end % for compatibility
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isfield(net, 'version'), version = net.version; else version = 1; end
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end

    % check compatibility
    if version == 1
        [DI, DIsub] = calcMvarDnnDI_(net, nodeControl, exControl, isFullNode);
        return;
    end

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc multivariate VAR DNN DI
    DI = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    nodeInputOrg = ones((nodeNum + exNum)*lags, 1);
    for i=1:nodeNum
        if isempty(net.nodeNetwork{i}), continue; end
        [~,idx] = find(control(i,:,:)==1);
        nodeInput = nodeInputOrg(idx,1);

        % predict 
        DIsub(i,1) = predict(net.nodeNetwork{i}, nodeInput, 'ExecutionEnvironment', 'cpu');

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            nodeInput = nodeInputOrg;
            for p=1:lags, nodeInput(j+(nodeNum + exNum)*(p-1),:) = 0; end
            nodeInput = nodeInput(idx,1);

            % predict 
            DIsub(i,j+1) = predict(net.nodeNetwork{i}, nodeInput, 'ExecutionEnvironment', 'cpu');
            DI(i,j) = abs(DIsub(i,1) - DIsub(i,j+1));
        end
    end
end

