%%
% Caluclate multivariate Granger Causality
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMultivariateGCI2(X, exSignal, nodeControl, exControl, lags, alpha, isFullNode)
    if nargin < 7
        isFullNode = 0;
    end
    if nargin < 6
        alpha = 0.05;
    end
    if nargin < 5
        lags = 3;
    end
    if nargin < 4
        exControl = [];
    end
    if nargin < 3
        nodeControl = [];
    end
    if nargin < 2
        exSignal = [];
    end
    nodeNum = size(X,1);
    nodeInNum = nodeNum + size(exSignal,1);
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % set node input
    if ~isempty(exSignal)
        X = [X; exSignal];
    end

    len = size(X,2);
    p = lags;
    Y = flipud(X.'); % need to flip signal

    % first, calculate multivariate autoregression without target
    Yjs = cell(1,nodeMax);
    nn1 = nodeMax-1;
    for j=1:nodeMax
        Ytj = zeros(len-p, p*nn1);
        Yj = [Y(:,1:j-1), Y(:,j+1:end)];
        for k=1:p
            Ytj(:,1+nn1*(k-1):nn1*k) = Yj(1+k:len-p+k,:);
        end
        Yjs{j} = Ytj;
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
        X2 = X;
        if ~isempty(nodeControl)
            filter = nodeControl(i,:).';
            X2(1:nodeNum,:) = X2(1:nodeNum,:) .* filter;
        end
        % input signal is time [1 ... last]
        % need to flip signal
        Xi = flipud(X2(i,:).');

        % multivariate autoregression
        Xt = Xi(1:len-p);
        Xti = zeros(len-p, p*nodeMax);
        for k=1:p
            Xti(:,1+nodeMax*(k-1):nodeMax*k) = Y(1+k:len-p+k,:);
        end
        % apply the regress function
        [b,bint,r] = regress(Xt,Xti);
        Vxt = var(r);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        T = len-p;
        RSS = r'*r;
        k = p * nodeMax;
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        for j=1:nodeMax
            if i==j, continue; end
            [b,bint,r] = regress(Xt,Yjs{j});
            Vyt = var(r);

            gcI(i,j) = log(Vyt / Vxt);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS1 = r'*r;
            k1 = p * nn1;
            AIC(i,j) = T*log(RSS1/T) + 2 * k1;
            BIC(i,j) = T*log(RSS1/T) + k1*log(T);

            % calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            %RSS1 = r'*r;  % p1 = p*nn1;
            RSS2 = RSS;   % p2 = p*nodeNum;
            F(i,j) = ((RSS1 - RSS2)/p) / (RSS2 / (len - k));
            P(i,j) = 1 - fcdf(F(i,j),p,(len-k));
            cvFd(i,j) = finv(1-alpha,p,(len-k));
            h(i,j) = F(i,j) > cvFd(i,j);
        end
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl(exControl==0) = nan;
        gcI(:,nodeNum+1:end) = gcI(:,nodeNum+1:end) .* exControl;
        F(:,nodeNum+1:end) = F(:,nodeNum+1:end) .* exControl;
        P(:,nodeNum+1:end) = P(:,nodeNum+1:end) .* exControl;
        cvFd(:,nodeNum+1:end) = cvFd(:,nodeNum+1:end) .* exControl;
        h(:,nodeNum+1:end) = h(:,nodeNum+1:end) .* exControl;
        AIC(:,nodeNum+1:end) = AIC(:,nodeNum+1:end) .* exControl;
        BIC(:,nodeNum+1:end) = BIC(:,nodeNum+1:end) .* exControl;
    end
end

