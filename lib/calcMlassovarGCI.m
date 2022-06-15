%%
% Calculate mLassoVAR (multivaliate Lasso Vector Auto-Regression) Granger Causality
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mLassoVAR network
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMlassovarGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    inputNum = nodeNum + net.exNum;
    lambda = net.lambda;
    elaAlpha = net.elaAlpha;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
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
        [~,idx] = find(control(i,:,:)==1);
        
        % vector auto-regression (VAR)
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);
        mc = length(find(abs(net.bvec{i})>0));

        % var of residuals of node
        r = Xt - (Xti*net.bvec{i} + net.info{i}.Intercept);
        Vxt = var(r,1);

        % AIC and BIC of this node (assuming residuals are gausiann distribution)
        T = sigLen-lags;
        RSS = r'*r;
        k = mc+1;
        nodeAIC(i) = T*log(RSS/T) + 2 * k;
        nodeBIC(i) = T*log(RSS/T) + k*log(T);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % PLS vector auto-regression (VAR)
            control2 = control;
            control2(i,j,:) = 0;
            [~,idx2] = find(control2(i,:,:)==1);
            Xtj = Yj(:,idx2);

            % apply the Lasso regress function
            [b,info] = lasso(Xtj,Xt,'Lambda',lambda,'Alpha',elaAlpha); % including Intercept
            r = Xt - (Xtj*b + info.Intercept);
            Vyt = var(r,1);

            gcI(i,j) = log(Vyt / Vxt);

            % AIC and BIC (assuming residuals are gausiann distribution)
            % BIC = n*ln(RSS/n)+k*ln(n)
            RSS1 = r'*r;
            mc = length(find(abs(b)>0));
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

