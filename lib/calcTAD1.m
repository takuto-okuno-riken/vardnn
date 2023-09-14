%%
% calculate TA-Δ1 (first-order temporal autocorrelation) (Shinn et al., 2022)
% returns TA-Δ1 (node x 1 values)
% input:
%  X         multivariate time series matrix (node x time series)

function TA = calcTAD1(X)
    R = corr(X(:,1:end-1)',X(:,2:end)');
    TA = diag(R);
end
