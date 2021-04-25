%%
% Caluclate mPLSVAR (multivaliate PLS Vector Auto-Regression) Granger Causality
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPLSVAR network
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMplsvarGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    nodeInNum = nodeNum + net.exNum;
    p = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end

    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-p, p*nodeInNum);
    for k=1:p
        Yj(:,1+nodeInNum*(k-1):nodeInNum*k) = Y(1+k:sigLen-p+k,:);
    end

    nodeAIC = zeros(nodeNum,1);
    nodeBIC = zeros(nodeNum,1);
    gcI = nan(nodeNum,nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeInNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:p
            idx = [idx, nodeIdx+nodeInNum*(k-1), exIdx+nodeInNum*(k-1)];
        end
        idxList = [nodeIdx, exIdx];
        nlen = length(idxList);
        Xt = Y(1:sigLen-p,i);
        Xti = Yj(:,idx);
        mc = net.ncomp;

        % var of residuals of node
%        r1 = Xt - [ones(size(Xti,1),1), Xti] * net.bvec{i}; % r1 == r
        r = net.stats{i}.Yresiduals;
        Vxt = var(r);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        T = sigLen-p;
        RSS = r'*r;
        k = mc+1;
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            Xtj = Xti;
            bIdx = find(idxList==j);
            for k=1:p
                Xtj(:,bIdx+nlen*(k-1)) = 0;
            end
            r = Xt - [ones(size(Xtj,1),1), Xtj] * net.bvec{i};
            Vyt = var(r);

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
            F(i,j) = ((RSS1 - RSS2)/p) / (RSS2 / (sigLen - k));
            P(i,j) = 1 - fcdf(F(i,j),p,(sigLen-k));
            cvFd(i,j) = finv(1-alpha,p,(sigLen-k));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
end

