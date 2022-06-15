%%
% Calculate time shifted Correlation (Functional Connectivity)
% cross-correlation type signal processing
% returns Functional Connectivity (FC)
% input:
%  X      multivariate time series matrix (node x time series)
%  p      number of lags for shifting

function [FC] = calcTimeShiftedCorrelation_(X, p)
    [n,m] = size(X);
    if p>0
        A = zeros(n,n,p);
        X1 = X(:,1:m-p);
        for i=1:p
            X2 = X(:,1+i:m-(p-i));
            A(:,:,i) = corr(X2.', X1.');
        end
        if p>1
            B = abs(A);
            [m,mI] = max(B,[],3); % get characterized index of each n x n mat
            C = repmat(1:n,n,1);
            D = repmat((0:n:(n-1)*n).',1,n);
            idx = (C+D).'+((mI-1)*(n*n)); % get index of characterized item of B
            Z = zeros(n,n,p);
            Z(idx(:)) = 1;
            A = A .* Z; % mask
            FC = sum(A,3);
        else
            FC = A;
        end
    else
        FC = corr(X.', X.');
    end
end
