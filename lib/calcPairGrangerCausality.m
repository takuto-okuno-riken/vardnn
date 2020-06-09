%%
% Caluclate pairwised Granger Causality
% returns Granger causality index (gcI)
% input:
%  X      time series vector (1 x time series)
%  Y      time series vector (1 x time series)
%  p      number of lags for autoregression

function gcI = calcPairGrangerCausality(X, Y, p)
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
    [b,bint,r] = regress(Xt,Xti);
    Vxt = var(r);

    % autoregression plus other regression
    Yt = X(1:n-p);
    Yti = zeros(n-p, p*2);
    for i=1:p
        Yti(:,i) = X(i+1:n-p+i);
        Yti(:,p+i) = Y(i+1:n-p+i);
    end
    [b,bint,r] = regress(Yt,Yti);
    Vyt = var(r);
    
    gcI = log(Vxt / Vyt);
end