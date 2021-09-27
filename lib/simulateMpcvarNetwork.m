%%
% multi-step ahead forecasting (recursive) node signals by PCVAR (Principal Component Vector Auto-Regression) and exogenous input
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPCVAR network

function [S, time] = simulateMpcvarNetwork(X, exSignal, nodeControl, exControl, net)
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    lags = net.lags;

    % set node input
    S = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    
    idxs = {};
    for i=1:nodeNum
        for k=1:lags
            [~,idxs{i,k}] = find(control(i,:,k)==1);
        end
    end

    disp('start multi-step ahead forecasting (recursive) by mPCVAR network');
    ticH = tic;
    maxComp = net.maxComp;
    nmu = net.mu;
    coeff = net.coeff;
    bvec = net.bvec;
    for t=lags:sigLen-1
        if mod(t,10)==0, disp(['step : ' num2str(t)]); end
        for i=1:nodeNum
            S2 = [];
            for k=1:lags
                S2 = [S2; S(idxs{i,k},t-(k-1))];
            end

            % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
            mc = maxComp{i};
            mu = nmu{i};

            % predict next time step
            score = (S2.' - mu) / coeff{i}.';
            subScore = [score(:,1:mc), 1];
            S(i,t+1) = subScore * bvec{i};
        end
        % fixed over shoot values
        idx = find(S(:,t+1) > 1.2);
        S(idx,t+1) = 1.2;
        idx = find(S(:,t+1) < -0.2);
        S(idx,t+1) = -0.2;
    end
    S = S(1:nodeNum,:);
    time = toc(ticH);
    disp(['finish multi-step ahead forecasting (recursive) by mPCVAR network! t = ' num2str(time) 's']);
end
