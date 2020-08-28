%%
% Simulate node signals by traind DLCM and exogenous input
% input:
%  X            multivariate time series matrix (node x time series)
%  inSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network

function [S, time] = simulateDlcmNetwork(X, inSignal, nodeControl, inControl, netDLCM)
    nodeNum = size(X,1);
    sigLen = size(X,2); % TODO:

    % set first signal
    S =  X(:,1:sigLen);

    disp('start simulation whole DLCM network');
    ticH = tic;
    for t=1:sigLen-1
        if mod(t,10)==0, disp(['step : ' num2str(t)]); end
        if isempty(inSignal)
            nodeInputOrg = S(:,t);
        else
            nodeInputOrg = [S(:,t); inSignal(:,t)];
        end
        for i=1:nodeNum
            nodeInput = nodeInputOrg;
            if ~isempty(nodeControl)
                filter = nodeControl(i,:).';
                nodeInput(1:nodeNum,1) = nodeInput(1:nodeNum,1) .* filter;
            end
            if ~isempty(inControl)
                filter = inControl(i,:).';
                nodeInput(nodeNum+1:end,1) = nodeInput(nodeNum+1:end,1) .* filter;
            end
            % predict next time step
            S(i,t+1) = predict(netDLCM.nodeNetwork{i}, nodeInput);
        end
        % fixed over shoot values
        idx = find(S(:,t+1) > 1.2);
        S(idx,t+1) = 1.2;
        idx = find(S(:,t+1) < -0.2);
        S(idx,t+1) = -0.2;
    end
    time = toc(ticH);
    disp(['finish simulation whole DLCM network! t = ' num2str(time) 's']);
end
