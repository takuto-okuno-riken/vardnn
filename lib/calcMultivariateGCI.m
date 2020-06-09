%%
% Caluclate multivariate Granger Causality
% returns Granger causality index matrix (gcI)
% VAR coefficients (A), VAR residuals (E) and VAR residuals covariance matrix (vE)
% input:
%  X      multivariate time series matrix (node x time series)
%  lags   number of lags for autoregression

function [gcI, A, E, vE] = calcMultivariateGCI(X, lags)
    n = size(X,1);
    gcI = nan(n,n);

    % full regression
    [A,E,vE]   = calcVARcoeff(X,lags);
    lvE = log(diag(vE)); % log of variances of residuals

    for j = 1:n
        Y = X;
        Y(j,:) = [];

        % reduced regression
        [~,~,rvE] = calcVARcoeff(Y,lags);
        lrvE = log(diag(rvE)); % log of variances of residuals

        G = nan(n,1);
        if j>1, G(1:j-1) = lrvE(1:j-1); end
        if j<n, G(j+1:n) = lrvE(j:n-1); end
        gcI(:,j) = G - lvE; % log rvE/vE
    end
end

%%
% Returns VAR coefficients (A), VAR residuals (E)
% and VAR residuals covariance matrix (vE)
% input:
%  X     multivariate time series matrix (node x time series)
%  p     number of lags for autoregression
function [A, E, vE] = calcVARcoeff(X, p)
    [n,m] = size(X);

    U = ones(1,m);
    X = X-mean(X,2)*U; % x - E(x) of time seriase
    rm = (m-p);

    % lags
    X0 = X(:,p+1:m); % unlagged observations
    Xl = zeros(n*p,rm);
    for k = 1:p
        s = n*(k-1);
        Xl(s+1:s+n,:) = X(:,p+1-k:m-k); % k-lagged observations
    end

    A  = X0/Xl; % OLS using QR decomposition (might include inf or nan)
    E  = X0-A*Xl;        % residuals
    vE = (E*E')/(rm-1);  % residuals covariance matrix

    A = reshape(A,n,n,p); % A(:,:,k) is the k-lag coefficients matrix
end
