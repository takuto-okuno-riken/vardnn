%%
% Calculate Pairwise VAR DNN Granger Causality
% returns Pairwise VAR DNN Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained pairwise VAR DNN network structure
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC] = calcPvarDnnGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end

    nodeNum = net.nodeNum;
    sigLen = size(X,2);
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + net.exNum; end

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    Y = flipud(Y.'); % need to flip signal

    % calc Pairwised VAR DNN granger causality
    gcI = nan(nodeNum, nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        [~,idx] = find(control(i,i,:)==1);
        Xt = Y(1:sigLen-lags,i);
        Yi = zeros(sigLen-lags, lags);
        for k=1:lags, Yi(:,k) = Y(1+k:sigLen-lags+k,i); end
        Xti = Yi(:,idx);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            [~,idx] = find(control(i,j,:)==1);
            Yj = zeros(sigLen-lags, lags);
            for k=1:lags, Yj(:,k) = Y(1+k:sigLen-lags+k,j); end
            Xtj = Yj(:,idx);

            % predict 
            Si = predict(net.nodeNetwork{i,j}, [Xti,Xtj].');
            err = Si - Xt.';
            VarEi = var(err,1);

            % AIC and BIC of this node (assuming residuals are gausiann distribution)
            T = sigLen-1;
            RSS = err*err';
            k = 2*lags + 1; % input + bias

            % imparement node signals
            Yj(:) = 0;
            Xtj = Yj(:,idx);

            % predict 
            Sj = predict(net.nodeNetwork{i,j}, [Xti,Xtj].');
            err = Sj - Xt.';
            VarEj = var(err,1);
            gcI(i,j) = log(VarEj / VarEi);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS1 = err*err';
            k1 = 2*lags + 1; % input + bias
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

