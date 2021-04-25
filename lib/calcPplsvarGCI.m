%%
% Caluclate pairwise PLSVAR Granger Causality
% returns pairwise PLSVAR Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained pairwise PLS VAR network structure
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC] = calcPplsvarGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end

    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;
    sigLen = size(X,2);
    p = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % set node input
    Y = [X; exSignal];

    % calc Pairwised PLS VAR granger causality
    gcI = nan(nodeNum, nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            Y1 = flipud(Y(i,:));
            Y2 = flipud(Y(j,:));

            % autoregression plus other regression
            Yt = Y2(1:sigLen-p).'; % TODO: X1 & X2 opposite ??
            Yti = zeros(sigLen-p, p*2);
            for k=1:p
                Yti(:,k) = Y2(k+1:sigLen-p+k);
%                Yti(:,p+k) = Y1(k+1:sigLen-p+k);
            end

            % var of residuals (full)
            r = net.stats{i,j}.Yresiduals;
            Vxt = var(r);

            % AIC and BIC of this node (assuming residuals are gausiann distribution)
            mc = net.ncomp;
            T = sigLen-p;
            RSS = r'*r;
            k = mc+1;

            % var of residuals (reduced)
            r = Yt - [ones(size(Yti,1),1), Yti] * net.bvec{i,j};
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

