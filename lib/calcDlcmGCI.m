%%
% Caluclate DLCM Granger Causality
% returns DLCM Granger causality index matrix (gcI) and AIC, BIC (of node) vector
% input:
%  X            multivariate time series matrix (node x time series)
%  inSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network

function [gcI, nodeAIC, nodeBIC] = calcDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM)
    nodeNum = size(X,1);
    sigLen = size(X,2);
    
    % set node input
    if isempty(inSignal)
        nodeInputOrg = X(:,1:sigLen-1);
    else
        nodeInputOrg = [X(:,1:sigLen-1); inSignal(:,1:sigLen-1)];
    end

    % calc DLCM granger causality
    nodeAIC = zeros(nodeNum,1);
    nodeBIC = zeros(nodeNum,1);
    gcI = nan(nodeNum, nodeNum);
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
        T = sigLen;
        RSS = err*err';
        k = nodeNum + size(inSignal, 1);
        for j=2:2:length(netDLCM.nodeNetwork{i, 1}.Layers)
            k = k + length(netDLCM.nodeNetwork{i, 1}.Layers(j, 1).Bias);
        end
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        % imparement node signals
        for j=1:nodeNum
            if i==j, continue; end
            impInput = nodeInput;
            impInput(j,:) = 0;
            % predict 
            Sj = predict(netDLCM.nodeNetwork{i}, impInput);
            err = Sj - nodeTeach;
            VarEj = var(err);
            gcI(i,j) = log(VarEj / VarEi);
        end
    end
end

