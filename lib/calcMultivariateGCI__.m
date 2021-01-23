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

function [gcI, A, E, vE] = calcMultivariateGCI__(X, exSignal, nodeControl, exControl, lags, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    nodeMax = nodeNum + size(exSignal,1);
    
    % set node input
    if ~isempty(exSignal)
        X = [X; exSignal];
    end

    gcI = nan(nodeMax,nodeMax);

    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeNum+size(exSignal,1)];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        % full regression
        Y = X([nodeIdx, exIdx],:);
        [A,E,vE]   = calcVARcoeff(Y,lags);
        rsd = log(diag(vE)); % log of variances of residuals
        lvE = nan(nodeMax,1);
        lvE([nodeIdx, exIdx]) = rsd;

        % reduced regression
        [~,idx] = find(nodeIdx==i); nodeIdx(idx) = [];
        Y = X([nodeIdx, exIdx],:);
        [~,~,rvE] = calcVARcoeff(Y,lags);
        lrvE = log(diag(rvE)); % log of variances of residuals

        G = nan(nodeMax,1);
        G([nodeIdx, exIdx]) = lrvE;
        gcI(:,i) = G - lvE; % log rvE/vE
    end
    gcI(1:nodeNum,nodeNum+1:end) = gcI(nodeNum+1:end, 1:nodeNum)';

    % output control
    if ~isempty(exSignal)
        gcI = gcI(1:nodeNum,:);
    end
    if isFullNode == 0
        gcI = gcI(:,1:nodeNum);
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
    X = X-mean(X,2)*U; % x - E(x) of time series
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
