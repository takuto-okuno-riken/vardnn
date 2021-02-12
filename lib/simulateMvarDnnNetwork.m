%%
% Simulate node signals by traind multivariate VAR DNN and exogenous input
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained mVAR DNN network

function [S, time] = simulateMvarDnnNetwork(X, exSignal, nodeControl, exControl, net)
    nodeNum = size(X,1);
    sigLen = size(X,2); % TODO:
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end

    % set first signal
    S = X;

    disp('start simulation whole mVAR DNN (DLCM) network');
    ticH = tic;

    for t=lags:sigLen-1
        if mod(t,10)==0, disp(['step : ' num2str(t)]); end
        nodeInputOrg = [];
        for i=1:lags, nodeInputOrg = [nodeInputOrg; S(:,t-lags+i)]; end
        if ~isempty(exSignal)
            for i=1:lags, nodeInputOrg = [nodeInputOrg; exSignal(:,t-lags+i)]; end
        end
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
            % predict next time step
            S(i,t+1) = predict(net.nodeNetwork{i}, nodeInput);
        end
        % fixed over shoot values
        idx = find(S(:,t+1) > 1.2);
        S(idx,t+1) = 1.2;
        idx = find(S(:,t+1) < -0.2);
        S(idx,t+1) = -0.2;
    end
    time = toc(ticH);
    disp(['finish simulation whole mVAR DNN (DLCM) network! t = ' num2str(time) 's']);
end
