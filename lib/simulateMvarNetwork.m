%%
% Simulate node signals by MVAR (multivaliate vector auto-regression) and exogenous input
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          MVAR network

function [S, time] = simulateMvarNetwork(X, exSignal, nodeControl, exControl, net)
    nodeNum = size(X,1);
    sigLen = size(X,2); % TODO:
    p = net.lags;

    % set node input
    S = [X; exSignal];

    disp('start simulation whole MVAR network');
    ticH = tic;
    for t=p:sigLen-1
        if mod(t,10)==0, disp(['step : ' num2str(t)]); end
        for i=1:nodeNum
            nodeIdx = [1:nodeNum];
            if ~isempty(nodeControl)
                [~,nodeIdx] = find(nodeControl(i,:)==1);
            end
            exIdx = [nodeNum+1:nodeNum+net.exNum];
            if ~isempty(exControl)
                [~,exIdx] = find(exControl(i,:)==1);
                exIdx = exIdx + nodeNum;
            end
            S2 = [];
            for k=1:p
                S2 = [S2; S([nodeIdx, exIdx],t-(k-1))];
            end
            S2 = [S2; 1]; % might not be good to add bias

            % predict next time step
            S(i,t+1) = S2.' * net.bvec{i};
        end
        % fixed over shoot values
        idx = find(S(:,t+1) > 1.2);
        S(idx,t+1) = 1.2;
        idx = find(S(:,t+1) < -0.2);
        S(idx,t+1) = -0.2;
    end
    S = S(1:nodeNum,:);
    time = toc(ticH);
    disp(['finish simulation whole MVAR network! t = ' num2str(time) 's']);
end
