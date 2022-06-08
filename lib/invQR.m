%%
% matrix inversion based on QR decomposition (accurate version)
% X should be Square matrix

function Xi = invQR(X)
    [sz1,sz2] = size(X);
    [Q,R,perm] = qr(X,0);
    p = sum(abs(diag(R)) > max(sz1,sz2)*eps(R(1))); % 2 steps differ than zero
    if p < sz2
        R = R(1:p,1:p);
        Q = Q(:,1:p);
        perm = perm(1:p);
    end
    Ci = inv(R) * Q';
    Xi = zeros(sz2,sz1);
    Xi(perm,:) = Ci;
end
