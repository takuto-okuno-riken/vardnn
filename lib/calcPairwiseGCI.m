%%
% Caluclate pairwise Granger Causality
% returns Granger causality index matrix (gcI)
% input:
%  X      multivariate time series matrix (node x time series)
%  lags   number of lags for autoregression

function gcI = calcPairwiseGCI(X, lags)
    n = size(X,1);
    gcI = nan(n,n);
    for i=1:n
        for j=1:n
            if i==j, continue; end
            gcI(i,j) = calcPairGrangerCausality(X(i,:), X(j,:), lags);
        end
    end
end

