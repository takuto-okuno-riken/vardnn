%%
% Caluclate DLCM Granger Causality
% returns DLCM Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  inSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI_(X, inSignal, nodeControl, inControl, netDLCM, alpha, isFullNode)
    if nargin < 7
        isFullNode = 0;
    end
    if nargin < 6
        alpha = 0.05;
    end
    nodeNum = size(X,1);
    nodeInNum = nodeNum + size(inSignal,1);
    sigLen = size(X,2);
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % set node input
    if isempty(inSignal)
        nodeInputOrg = X(:,1:sigLen-1);
    else
        nodeInputOrg = [X(:,1:sigLen-1); inSignal(:,1:sigLen-1)];
    end

    % calc DLCM granger causality
    nodeAIC = zeros(nodeNum,1);
    nodeBIC = zeros(nodeNum,1);
    gcI = nan(nodeNum, nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);
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
        nodeTeach = X(i,2:end);
        % predict 
        Si = predict(netDLCM.nodeNetwork{i}, nodeInput);
        err = Si - nodeTeach;
        VarEi = var(err);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        T = sigLen-1;
        RSS = err*err';
        k = nodeNum + size(inSignal, 1) + 1; % input + bias
        %for j=2:2:length(netDLCM.nodeNetwork{i, 1}.Layers)
        %    k = k + length(netDLCM.nodeNetwork{i, 1}.Layers(j, 1).Bias);   % added hidden neuron number
        %end
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            impInput = nodeInput;
            impInput(j,:) = 0;

            % predict 
            Sj = predict(netDLCM.nodeNetwork{i}, impInput);
            err = Sj - nodeTeach;
            VarEj = var(err);
            gcI(i,j) = log(VarEj / VarEi);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS1 = err*err';
            k1 = nodeNum - 1 + size(inSignal, 1) + 1;
            AIC(i,j) = T*log(RSS1/T) + 2 * k1;
            BIC(i,j) = T*log(RSS1/T) + k1*log(T);

            % calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            %RSS1 = err*err';  % p1 = nodeNum - 1 + size(inSignal, 1) + 1;
            RSS2 = RSS;       % p2 = k
            F(i,j) = ((RSS1 - RSS2)/1) / (RSS2 / (sigLen - k));
            P(i,j) = 1 - fcdf(F(i,j),1,(sigLen - k));
            cvFd(i,j) = finv(1-alpha,1,(sigLen - k));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
end

