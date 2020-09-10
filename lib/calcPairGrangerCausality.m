%%
% Caluclate pairwised Granger Causality
% returns Granger causality index (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F) and the critical value from the F-distribution (cvFd)
% input:
%  X      time series vector (1 x time series)
%  Y      time series vector (1 x time series)
%  p      number of lags for autoregression
%  alpha  the significance level of F-statistic (default:0.05)

function [gcI, h, P, F, cvFd] = calcPairGrangerCausality(X, Y, p, alpha)
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
    Xti = zeros(n-p, p);
    for i=1:p
        Xti(:,i) = X(i+1:n-p+i);
    end
    % apply the regress function
    [b,bint,Xr] = regress(Xt,Xti);
    Vxt = var(Xr);

    % autoregression plus other regression
    Yt = X(1:n-p);
    Yti = zeros(n-p, p*2);
    for i=1:p
        Yti(:,i) = X(i+1:n-p+i);
        Yti(:,p+i) = Y(i+1:n-p+i);
    end
    [b,bint,Yr] = regress(Yt,Yti);
    Vyt = var(Yr);
    
    gcI = log(Vxt / Vyt);

    % calc F-statistic
    % https://en.wikipedia.org/wiki/F-test
    % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
    RSS1 = Xr'*Xr;
    RSS2 = Yr'*Yr;
    F = ((RSS1 - RSS2)/p) / (RSS2 / (n - 2*p));
    P = 1 - fcdf(F,p,(n-2*p));
    cvFd = finv(1-alpha,p,(n-2*p));
    h = F > cvFd;
end