%%
% Calculate multivariate VAR DNN Granger causality
% returns multivariate VAR DNN Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained multivariate VAR DNN network
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMvarDnnGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isfield(net, 'version'), version = net.version; else version = 1; end

    % check compatibility
    if version == 1
        [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMvarDnnGCI_(X, exSignal, nodeControl, exControl, net, alpha, isFullNode);
        return;
    end

    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    sigLen = size(X,2);
    inputNum = nodeNum + exNum;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    Y = flipud(Y.'); % need to flip signal

    % set training data
    Yj = zeros(sigLen-lags, lags*inputNum);
    for p=1:lags
        Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sigLen-lags+p,:);
    end

    % calc mVAR DNN granger causality
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
        if isempty(net.nodeNetwork{i}), continue; end
        [~,idx] = find(control(i,:,:)==1);
        
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);

        % predict 
        Si = predict(net.nodeNetwork{i}, Xti.', 'ExecutionEnvironment', 'cpu');
        err = Si - Xt.';
        VarEi = var(err,1);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        T = sigLen-1;
        RSS = err*err';
        k = nodeNum + size(exSignal, 1) + 1; % input + bias
        %for j=2:2:length(net.nodeNetwork{i, 1}.Layers)
        %    k = k + length(net.nodeNetwork{i, 1}.Layers(j, 1).Bias);   % added hidden neuron number
        %end
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            
            Xtj = Yj;
            for p=1:lags, Xtj(:,j+(nodeNum + exNum)*(p-1)) = 0; end
            Xtj = Xtj(:,idx);

            % predict 
            Sj = predict(net.nodeNetwork{i}, Xtj.', 'ExecutionEnvironment', 'cpu');
            err = Sj - Xt.';
            VarEj = var(err,1);
            gcI(i,j) = log(VarEj / VarEi);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS1 = err*err';
            k1 = nodeNum - 1 + size(exSignal, 1) + 1;
            AIC(i,j) = T*log(RSS1/T) + 2 * k1;
            BIC(i,j) = T*log(RSS1/T) + k1*log(T);

            % calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            %RSS1 = err*err';  % p1 = nodeNum - 1 + size(exSignal, 1) + 1;
            RSS2 = RSS;       % p2 = k
            F(i,j) = ((RSS1 - RSS2)/1) / (RSS2 / (sigLen - k));
            P(i,j) = 1 - fcdf(F(i,j),1,(sigLen - k));
            cvFd(i,j) = finv(1-alpha,1,(sigLen - k));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
end

