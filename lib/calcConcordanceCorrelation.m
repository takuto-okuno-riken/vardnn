%%
% Calculate Concordance Correlation Coefficient
% L. I-K. Lin, Biometrics, Vol. 45, No. 1 (Mar., 1989), pp. 255-268
% returns CCC (CCC), P is corr result.
% input:
%  X            multivariate time series matrix (node x time series)

function [CCC, P] = calcConcordanceCorrelation(X)
    n = size(X,1);
    [CM, P] = corr(X.', X.');
    M = repmat(mean(X,2),1,n);  % mean to matrix
    S = repmat(std(X,1,2),1,n); % sigma to matrix
    V = S ./ S'; % scale shift
    U = (M - M') ./ sqrt(S .* S'); % location shift
    CB = ((V + 1./V + U.^2)/2).^-1;
    CCC = CM .* CB;
end
