%%
% Caluclate multivariate VAR LSTM Directional Influence matrix (DI) and impaired node signals (DIsub)
% returns mVAR LSTM DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained multivariate VAR LSTM network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = calcMvarDnnDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end

    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    if isfield(net, 'exNum'), exNum = net.exNum; else exNum = size(net.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2) - nodeNum; end % for compatibility
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end

    % calc multivariate VAR LSTM DI
    DI = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    for i=1:nodeNum
        nodeInput = ones(nodeNum + exNum, lags);
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', 1, lags);
            nodeInput(1:nodeNum,:) = nodeInput(1:nodeNum,:) .* filter;
        end
        if ~isempty(exControl)
            filter = repmat(exControl(i,:).', 1, lags);
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        % predict 
        DIsub(i,1)  = predict(net.nodeNetwork{i}, nodeInput, 'ExecutionEnvironment', 'cpu');

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            impInput = nodeInput;
            impInput(j,:) = 0;

            % predict 
            DIsub(i,j+1) = predict(net.nodeNetwork{i}, impInput, 'ExecutionEnvironment', 'cpu');
            DI(i,j) = abs(DIsub(i,1) - DIsub(i,j+1));
        end
    end
end

