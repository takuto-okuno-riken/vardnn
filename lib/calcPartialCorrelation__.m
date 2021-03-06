%%
% Caluclate Partial Correlation by Inverse Covariance
% returns Partial Correlation (PC) and p-values (P)
% input:
%  X       multivariate time series matrix (node x time series)
%  lambda  regularisation level

function [PC] = calcPartialCorrelation__(X)
    n = size(X,1);
    m = size(X,2);
    PC = nan(n,n);
    Y = 1:n;
    for i=1:n
        for j=1:n
            x = X(i,:).';
            y = X(j,:).';
            Yij = setdiff(setdiff(Y,i),j);
            z = X(Yij,:).';

            [b1,bint1,r1] = regress(x,[z, ones(m,1)]);
            [b2,bint2,r2] = regress(y,[z, ones(m,1)]);
            PC(i,j) = (r1.'*r2) / (sqrt(r1.'*r1)*sqrt(r2.'*r2));
        end
    end
end
