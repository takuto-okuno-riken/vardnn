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
    inputNum = nodeNum + net.exNum;
    sigLen = size(X,2);
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end
    
    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc Pairwise PLS VAR granger causality
    gcI = nan(nodeNum, nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);

    mc = net.ncomp;
    T = sigLen-lags;
    k1 = mc+1;
    for i=1:nodeNum
        [~,idx] = find(control(i,i,:)==1);
        Xt = Y(1:sigLen-lags,i);
        Yi = zeros(sigLen-lags, lags);
        for k=1:lags, Yi(:,k) = Y(1+k:sigLen-lags+k,i); end
        Xti = Yi(:,idx);

        % find component number
        ncomp = floor(size(Xti,2) / 2);
        if ncomp < 1, ncomp = 1; end
        if ncomp > 50, ncomp = 50; end

        % apply the PLS regress function
        [XL,YL,XS,YS,b,PCTVAR,MSE,stats] = plsregress(Xti,Xt,ncomp);
        r = stats.Yresiduals;
        Vxt = var(r,1);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        RSS1 = r'*r;

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % var of residuals (full)
            r = net.stats{i,j}.Yresiduals;
            Vyt = var(r,1);

            gcI(i,j) = log(Vxt / Vyt);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS2 = r'*r;
            AIC(i,j) = T*log(RSS2/T) + 2 * k1;
            BIC(i,j) = T*log(RSS2/T) + k1*log(T);

            % calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            F(i,j) = ((RSS1 - RSS2)/lags) / (RSS2 / (sigLen - k1));
            P(i,j) = 1 - fcdf(F(i,j),lags,(sigLen-k1));
            cvFd(i,j) = finv(1-alpha,lags,(sigLen-k1));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
end

