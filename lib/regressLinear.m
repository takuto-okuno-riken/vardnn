%%
% calculate linear regression.
% returns coefficient vector (b), residuals vector (r), T-value vector (T)
% and p-value vector (P)
% out
% input:
%  X       observations x regressors matrix
%  y       response vector

function [b, r, T, P, df, s, se] = regressLinear(y, X)
    X2 = (X'*X);
    % pinv() is slower than inv().
    % but it can avoid singular matrix and nearly singular matrix
    X2i = pinv(X2);
    b = X2i * X' * y;
    r  = y - X*b;
    df = size(X,1) - size(X,2); % degree of freedom
    s  = sqrt(sum(r.*r)/df);
    se = sqrt(diag(X2i)*s^2);
    T  = b./se;
    P  = (T>=0).*(1 - tcdf(T,df))*2 + (T<0).*(tcdf(T,df))*2;
end