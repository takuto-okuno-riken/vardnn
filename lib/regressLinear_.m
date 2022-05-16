%%
% calculate linear regression (faster version).
% returns coefficient vector (b), residuals vector (r), T-value vector (T)
% and p-value vector (P)
% out
% input:
%  X       observations x regressors matrix
%  y       response vector
%  X2i     inv(X'*X) matrix (optional)
%  dX2i    diag(X2i) (optional)

function [b, r, T, P, df, s, se] = regressLinear_(y, X, X2i, dX2i)
    if nargin < 4, dX2i = []; end
    if nargin < 3, X2i = []; end
    if isempty(X2i)
        X2 = (X'*X);
        % pinv() is slower than inv().
        % but it can avoid singular matrix and nearly singular matrix.
        % this calculation might take time in huge matrix.
        X2i = pinv(X2);
    end
    if isempty(dX2i)
        dX2i = diag(X2i);
    end
    
    % this () is important for speed, but it might affect T-value by floating point error
    % b and r are not so affected
    b  = X2i * (X' * y);
    r  = y - X*b;
    df = size(X,1) - size(X,2); % degree of freedom
    s  = sum(r.*r)/df;
    se = sqrt(dX2i*s);
    T  = b./se;
    dT = single(T); % to avoid 'half' error
    P  = (dT>=0).*(1 - tcdf(dT,df))*2 + (dT<0).*(tcdf(dT,df))*2;
end
