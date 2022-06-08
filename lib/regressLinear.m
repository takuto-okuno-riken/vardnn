%%
% calculate linear regression (faster version).
% returns coefficient vector (b), residuals vector (r), T-value vector (T)
% and p-value vector (P)
% input:
%  X       observations x regressors matrix
%  y       response vector
%  Q       Q of QR decomposition of X (optional)
%  R       R of QR decomposition of X (optional)
%  perm    perm of QR decomposition of X (optional)

function [b, r, T, P, df, s, se, rsq] = regressLinear(y, X, Q, R, perm, RiQ, dR2i)
    if nargin < 7, dR2i = []; end
    if nargin < 6, RiQ = []; end
    if nargin < 5, perm = []; end
    if nargin < 4, R = []; end
    if nargin < 3, Q = []; end

    % QR decomposition of X
    [sz1,sz2] = size(X);
    if isempty(R) && isempty(perm) && (isempty(RiQ) || isempty(dR2i))
        [Q,R,perm] = qr(X,0);
        if isempty(R)
            p = 0;
        else
            p = sum(abs(diag(R)) > max(sz1,sz2)*eps(R(1))); % 2 steps differ than zero
        end
        if p < sz2
            R = R(1:p,1:p);
            Q = Q(:,1:p);
            perm = perm(1:p);
        end
    end

    % output coefficient and residuals
    b = zeros(sz2,1);
    % this () is important for speed, but it might affect T-value by floating point error
    % b and r are not so affected
    if isempty(RiQ)
        b(perm) = R \ (Q'*y);
    else
        b(perm) = RiQ *y;
    end
    r  = y - X*b;

    % output T-values
    if nargout >= 3
        if isempty(dR2i)
            % QR decomposition is better accuracy than pinv() or just inv().
            R2i = invQR(R'*R);
            dR2i = diag(R2i);
        end
        dX2i = ones(sz2,1);
        dX2i(perm) = dR2i;
    
        df = sz1 - sz2; % degree of freedom
        rs = single(r); % 'half' may not enough accurate
        s  = sum(rs.*rs)/df;
        se = sqrt(dX2i*s);
        T  = b./se;
    end
    % output P-values
    if nargout >= 4
        dT = single(T); % to avoid 'half' error
        P  = (dT>=0).*(1 - tcdf(dT,df))*2 + (dT<0).*(tcdf(dT,df))*2;
    end
    % output R-squared
    if nargout >= 8
        SSres = sum(rs.*rs);
        ydf = y - mean(y);
        SStot = sum(ydf.*ydf);
        rsq = 1 - SSres / SStot;
    end
end

