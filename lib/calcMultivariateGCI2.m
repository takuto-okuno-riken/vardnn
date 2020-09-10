%%
% Caluclate multivariate Granger Causality
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node) vector
% input:
%  X      multivariate time series matrix (node x time series)
%  lags   number of lags for autoregression
%  alpha  the significance level of F-statistic (default:0.05)

function [gcI, h, P, F, cvFd, nodeAIC, nodeBIC] = calcMultivariateGCI2(X, lags, alpha)
    if nargin < 3
        alpha = 0.05;
    end
    nodeNum = size(X,1);
    n = size(X,2);
    p = lags;
    Y = flipud(X.'); % need to flip signal

    % first, calculate multivariate autoregression without target
    Yjs = cell(1,nodeNum);
    nn1 = nodeNum-1;
    for j=1:nodeNum
        Ytj = zeros(n-p, p*nn1);
        Yj = [Y(:,1:j-1), Y(:,j+1:end)];
        for k=1:p
            Ytj(:,1+nn1*(k-1):nn1*k) = Yj(1+k:n-p+k,:);
        end
        Yjs{j} = Ytj;
    end

    nodeAIC = zeros(nodeNum,1);
    nodeBIC = zeros(nodeNum,1);
    gcI = nan(nodeNum,nodeNum);
    h = nan(nodeNum,nodeNum);
    P = nan(nodeNum,nodeNum);
    F = nan(nodeNum,nodeNum);
    cvFd = nan(nodeNum,nodeNum);
    for i=1:nodeNum
        % input signal is time [1 ... last]
        % need to flip signal
        Xi = flipud(X(i,:).');

        % multivariate autoregression
        Xt = Xi(1:n-p);
        Xti = zeros(n-p, p*nodeNum);
        for k=1:p
            Xti(:,1+nodeNum*(k-1):nodeNum*k) = Y(1+k:n-p+k,:);
        end
        % apply the regress function
        [b,bint,r] = regress(Xt,Xti);
        Vxt = var(r);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        T = n-p;
        RSS = r'*r;
        k = p * nodeNum;
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        for j=1:nodeNum
            if i==j, continue; end
            [b,bint,r] = regress(Xt,Yjs{j});
            Vyt = var(r);

            gcI(i,j) = log(Vyt / Vxt);
            
            % calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            RSS1 = r'*r;  % p1 = p*nn1;
            RSS2 = RSS;   % p2 = p*nodeNum;
            F(i,j) = ((RSS1 - RSS2)/p) / (RSS2 / (n - (p*nodeNum)));
            P(i,j) = 1 - fcdf(F(i,j),p,(n-p*nodeNum));
            cvFd(i,j) = finv(1-alpha,p,(n-p*nodeNum));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
end

