%%
% Caluclate multivariate Granger Causality
% returns Granger causality index matrix (gcI)
% VAR coefficients (A), VAR residuals (E) and VAR residuals covariance matrix (vE)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, A, E, vE] = calcMultivariateGCI(X, exSignal, nodeControl, exControl, lags, isFullNode)
    if nargin < 6
        isFullNode = 0;
    end
    if nargin < 5
        lags = 3;
    end
    if nargin < 4
        exControl = [];
    end
    if nargin < 3
        nodeControl = [];
    end
    if nargin < 2
        exSignal = [];
    end
    nodeNum = size(X,1);
    nodeInNum = nodeNum + size(exSignal,1);
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % set node input
    if ~isempty(exSignal)
        X = [X; exSignal];
    end

    gcI = nan(nodeNum,nodeMax);

    % full regression
    [A,E,vE]   = calcVARcoeff(X,lags);
    lvE = log(diag(vE)); % log of variances of residuals

    for j = 1:nodeMax
        Y = X;
        Y(j,:) = [];

        % reduced regression
        [~,~,rvE] = calcVARcoeff(Y,lags);
        lrvE = log(diag(rvE)); % log of variances of residuals

        G = nan(nodeMax,1);
        if j>1, G(1:j-1) = lrvE(1:j-1); end
        if j<nodeMax, G(j+1:nodeMax) = lrvE(j:nodeMax-1); end
        gcI(:,j) = G - lvE(1:nodeMax); % log rvE/vE
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        gcI(:,1:nodeNum) = gcI(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        gcI(:,nodeNum+1:end) = gcI(:,nodeNum+1:end) .* exControl;
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
