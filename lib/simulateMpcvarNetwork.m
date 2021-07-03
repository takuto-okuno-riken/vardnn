%%
% Simulate node signals by PCVAR (Principal Component Vector Auto-Regression) and exogenous input
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPCVAR network

function [S, time] = simulateMpcvarNetwork(X, exSignal, nodeControl, exControl, net)
    nodeNum = size(X,1);
    sigLen = size(X,2); % TODO:
    p = net.lags;

    % set node input
    S = [X; exSignal];

    % set node control index
    idxs = {};
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
        idxs{i} = [nodeIdx, exIdx];
    end

    disp('start simulation whole mPCVAR network');
    ticH = tic;
    maxComp = net.maxComp;
    nmu = net.mu;
    coeff = net.coeff;
    bvec = net.bvec;
    for t=p:sigLen-1
        if mod(t,10)==0, disp(['step : ' num2str(t)]); end
        S3 = S(:,t+1);
        for i=1:nodeNum
%        parfor i=1:nodeNum
            S2 = [];
            for k=1:p
                S2 = [S2; S(idxs{i},t-(k-1))];
            end

            % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
            mc = maxComp{i};
            mu = nmu{i};

            % predict next time step
            score = (S2.' - mu) / coeff{i}.';
            subScore = [score(:,1:mc), 1];
            S3(i) = subScore * bvec{i};
        end
        S(:,t+1) = S3;
        % fixed over shoot values
        idx = find(S(:,t+1) > 1.2);
        S(idx,t+1) = 1.2;
        idx = find(S(:,t+1) < -0.2);
        S(idx,t+1) = -0.2;
    end
    S = S(1:nodeNum,:);
    time = toc(ticH);
    disp(['finish simulation whole mPCVAR network! t = ' num2str(time) 's']);
end
