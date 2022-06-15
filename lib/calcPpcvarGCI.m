%%
% Calculate pairwise PCVAR Granger Causality
% returns pairwise PCVAR Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained pairwise PCVAR network structure
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC] = calcPpcvarGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end

    nodeNum = net.nodeNum;
    sigLen = size(X,2);
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + net.exNum; end

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc Pairwised PCVAR granger causality
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
            if j<=nodeNum && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && exControl(i,j-nodeNum) == 0, continue; end

            % autoregression plus other regression
            [~,idx] = find(control(i,j,:)==1);
            Yj = zeros(sigLen-lags, lags);
            for k=1:lags, Yj(:,k) = Y(1+k:sigLen-lags+k,j); end
            Xtj = [Xti, Yj(:,idx)];

            % var of residuals (full)
            r = net.rvec{i,j};
            Vxt = var(r,1);

            % AIC and BIC of this node (assuming residuals are gausiann distribution)
            mc = net.maxComp{i,j};
            mu = net.mu{i,j};
            T = sigLen-lags;
            RSS = r'*r;
            k = mc+1;

            % var of residuals (reduced)
            scorej = (Xtj - mu) / net.coeff{i,j}.';
            pcXti = [scorej(:,1:mc), ones(sigLen-lags,1)]; % might not be good to add bias
            r = Xt - pcXti * net.bvec{i,j}; % calc residuals
            Vyt = var(r,1);

            gcI(i,j) = log(Vyt / Vxt);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS1 = r'*r;
            k1 = mc+1;
            AIC(i,j) = T*log(RSS1/T) + 2 * k1;
            BIC(i,j) = T*log(RSS1/T) + k1*log(T);

            % calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            %RSS1 = r'*r;  % p1 = p*nn1;
            RSS2 = RSS;   % p2 = p*nodeNum;
            F(i,j) = ((RSS1 - RSS2)/lags) / (RSS2 / (sigLen - k));
            P(i,j) = 1 - fcdf(F(i,j),lags,(sigLen-k));
            cvFd(i,j) = finv(1-alpha,lags,(sigLen-k));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
end

