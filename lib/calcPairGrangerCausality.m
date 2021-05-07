%%
% Caluclate pairwised Granger Causality
% returns Granger causality index (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F) and the critical value from the F-distribution (cvFd)
% and AIC, BIC
% input:
%  X      time series vector (1 x time series)
%  Y      time series vector (1 x time series)
%  p      number of lags for autoregression
%  alpha  the significance level of F-statistic (default:0.05)

function [gcI, h, P, F, cvFd, AIC, BIC] = calcPairGrangerCausality(X, Y, p, alpha)
    if nargin < 4
        alpha = 0.05;
    end
    % input signal is time [1 ... last]
    % need to inversion signal
    X = flipud(X(:));
    Y = flipud(Y(:));
    n = length(X);

    % autoregression
    Xt = X(1:n-p);
    Xti = ones(n-p, p+1);
    for i=1:p
        Xti(:,i+1) = X(i+1:n-p+i);
    end
    % apply the regress function
    [b,bint,Xr] = regress(Xt,Xti);
    Vxt = var(Xr,1);

    % autoregression plus other regression
    Yt = X(1:n-p);
    Yti = ones(n-p, p*2+1);
    for i=1:p
        Yti(:,i+1) = X(i+1:n-p+i);
        Yti(:,p+i+1) = Y(i+1:n-p+i);
    end
    [b,bint,Yr] = regress(Yt,Yti);
    Vyt = var(Yr,1);
    
    gcI = log(Vxt / Vyt);

    % AIC and BIC (assuming residuals are gausiann distribution)
    T = n-p;
    RSS = Yr'*Yr;
    k = p*2 + 1;
    AIC = T*log(RSS/T) + 2 * k;
    BIC = T*log(RSS/T) + k*log(T);

    % calc F-statistic
    % https://en.wikipedia.org/wiki/F-test
    % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
    RSS1 = Xr'*Xr;  % p1 = p + 1
    RSS2 = RSS;     % p2 = k
    F = ((RSS1 - RSS2)/p) / (RSS2 / (n - k));
    P = 1 - fcdf(F,p,(n-k));
    cvFd = finv(1-alpha,p,(n-k));
    h = F > cvFd;
end