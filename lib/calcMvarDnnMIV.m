%%
% Caluclate multivariate VAR DNN Mean Impact Value (MIV)
% returns multivariate VAR DNN Mean Impact Value matrix (MIV) and Mean Absolute Impact Value matrix (MAIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained multivariate VAR DNN network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [MIV, MAIV] = calcMvarDnnMIV(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end

    nodeNum = size(X,1);
    nodeInNum = nodeNum + size(exSignal,1);
    sigLen = size(X,2);
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end

    % set node input
    nodeInputOrg = [];
    for i=1:lags, nodeInputOrg = [nodeInputOrg; X(:,i:end-(lags-i+1))]; end
    for i=1:lags, nodeInputOrg = [nodeInputOrg; exSignal(:,i:end-(lags-i+1))]; end

    % calc mVAR DNN MIV
    MIV = nan(nodeNum, nodeMax);
    MAIV = nan(nodeNum, nodeMax);
    for i=1:nodeNum
        nodeInput = nodeInputOrg;
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', lags, size(nodeInput,2));
            nodeInput(1:nodeNum*lags,:) = nodeInput(1:nodeNum*lags,:) .* filter;
        end
        if ~isempty(exControl)
            filter = repmat(exControl(i,:).', lags, size(nodeInput,2));
            nodeInput(nodeNum*lags+1:end,:) = nodeInput(nodeNum*lags+1:end,:) .* filter;
        end
        nodeTeach = X(i,1+lags:end);

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            Xj1 = nodeInput;
            Xj2 = nodeInput;
            for p=1:lags
                Xj1(j+nodeNum*(p-1),:) = Xj1(j+nodeNum*(p-1),:) * 1.1; 
                Xj2(j+nodeNum*(p-1),:) = Xj2(j+nodeNum*(p-1),:) * 0.9; 
            end

            % predict 
            N1 = predict(net.nodeNetwork{i}, Xj1, 'ExecutionEnvironment', 'cpu');
            N2 = predict(net.nodeNetwork{i}, Xj2, 'ExecutionEnvironment', 'cpu');
            MIV(i,j) = mean(N1 - N2);
            MAIV(i,j) = mean(abs(N1 - N2));
        end
    end
end

