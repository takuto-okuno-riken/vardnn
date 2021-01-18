%%
% Caluclate Partial Correlation by Inverse Covariance
% returns Partial Correlation (PC) and p-values (P)
% input:
%  X       multivariate time series matrix (node x time series)
%  lambda  regularisation level

function [iCov] = calcPartialCorrelation_(X,lambda)
    if nargin < 2
        lambda = 0;
    end
    sigma = cov(X.');
    n = size(sigma,1);
%    B = L1precisionBCD(sigma,lambda);
    B = inv(sigma);
%    B = sigma\eye(n);
    iCov = nan(n,n);
    for i=1:n
        for j=1:n
            iCov(i,j) = -B(i,j) / sqrt(B(i,i)*B(j,j));
        end
    end
end
