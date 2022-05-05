%%
% calculate linear regression.
% returns coefficient vector (b), residuals vector (r), T-value vector (T)
% and p-value vector (P)
% out
% input:
%  X       observations x regressors matrix
%  y       response vector

function [b, r, T, P, df, s, se] = regressLinear(y, X)
    b = (X'*X)\X'*y;

    df = size(X,1) - size(X,2); % degree of freedom
    r  = y - X*b;
    s  = sqrt(sum(r.*r)/df);
    se = sqrt(diag(inv(X'*X))*s^2);
    T  = b./se;
    P  = (T>=0).*(1 - tcdf(T,df))*2 + (T<0).*(tcdf(T,df))*2;
end