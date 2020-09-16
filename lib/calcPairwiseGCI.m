%%
% Caluclate pairwise Granger Causality
% returns Granger causality index (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F) and the critical value from the F-distribution (cvFd)
% and AIC, BIC
% input:
%  X      multivariate time series matrix (node x time series)
%  lags   number of lags for autoregression
%  alpha  the significance level of F-statistic (default:0.05)

function [gcI, h, P, F, cvFd, AIC, BIC] = calcPairwiseGCI(X, lags, alpha)
    if nargin < 3
        alpha = 0.05;
    end
    n = size(X,1);
    gcI = nan(n,n);
    h = nan(n,n);
    P = nan(n,n);
    F = nan(n,n);
    cvFd = nan(n,n);
    AIC = nan(n,n);
    BIC = nan(n,n);
    for i=1:n
        for j=1:n
            if i==j, continue; end
            [gcI(i,j), h(i,j), P(i,j), F(i,j), cvFd(i,j), AIC(i,j), BIC(i,j)] = calcPairGrangerCausality(X(i,:), X(j,:), lags, alpha);
        end
    end
end

