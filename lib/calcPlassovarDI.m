%%
% Caluclate pLassoVAR (pairwise Lasso Vector Auto-Regression) DI
% returns pLassoVAR DI matrix (DI) and impaired node signals (DIsub)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  lambda       lambda for the Lasso (default:0.01)
%  elaAlpha     Elastic Net Alpha for the Lasso (default:1)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub, coeff] = calcPlassovarDI(X, exSignal, nodeControl, exControl, lags, lambda, elaAlpha, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, elaAlpha = 1; end
    if nargin < 6, lambda = 0.01; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    p = lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + size(exSignal,1); end
    
    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal
    
    % calc Pairwise LassoVAR
    DI = nan(nodeNum, nodeMax);
    coeff = nan(nodeNum, nodeMax);
    DIsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            % autoregression plus other regression
            Yt = Y(1:(sigLen-p),i); % target
            Yti = ones(sigLen-p, p*2);
            for k=1:p
                Yti(:,k) = Y(k+1:sigLen-p+k,i); % target
                Yti(:,p+k) = Y(k+1:sigLen-p+k,j); % source
            end
            [b,~] = lasso(Yti,Yt,'Lambda',lambda,'Alpha',elaAlpha);  % including Intercept
            DIsub(i,j,1) = sum(b);
            DIsub(i,j,2) = sum(b(1:p));
            coeff(i,j) = DIsub(i,j,1)-DIsub(i,j,2); % actually this is sum of b(p+2:end)
            DI(i,j) = abs(coeff(i,j));
        end
    end
end

