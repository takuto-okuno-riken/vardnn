%%
% preparation for linear regression (faster version).
% returns Q, R, perm, RiQ and dR2i
% input:
%  X       observations x regressors matrix

function [Q, R, perm, RiQ, dR2i] = regressPrepare(X)
    % QR decomposition of X
    if isa(X,'half')
        [Q, R, perm] = qr(single(X), 0);
    else
        [Q, R, perm] = qr(X, 0);
    end
    [sz1,sz2] = size(X);
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
    if nargout >= 4
        if size(R,1) == size(R,2)
            RiQ = inv(R) * Q';
        else
            RiQ = [];
        end
    end
    if nargout >= 5
        R2i = invQR(R'*R);
        dR2i = diag(R2i);
    end
end
