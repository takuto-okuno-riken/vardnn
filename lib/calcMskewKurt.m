%%
% calculate multivariate skewness and kurtosis based on Mardia (1970).
% returns multivariate skewness and kurtosis
% input:
%  X         multivariate time series matrix (node x time series)
%  C         covariance matrix (optional)

function [mskew, mkurto, D] = calcMskewKurt(X, C)
    if nargin < 2, C = []; end
    n = size(X,2);
    X = X' - mean(X',1);
    if isa(X,'half')
        X = single(X);
    end
    if isempty(C)
        C = cov(X,1);
    end
    if det(C) == 0
        iC = invQR(C);
    else
        iC = inv(C);
    end
    D = X * iC * X';
    
    mskew = sum(D.^3,'all')/(n*n); % multivariate skewness
    dD = diag(D);
    mkurto = dD'*dD/n;  % multivariate kurtosis
end
