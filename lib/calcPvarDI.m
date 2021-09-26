%%
% Caluclate pVAR (pairwise Vector Auto-Regression) DI
% returns pVAR DI matrix (DI) and impaired node signals (DIsub)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub, coeff] = calcPvarDI(X, exSignal, nodeControl, exControl, lags, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    p = lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end
    
    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc Pairwise VAR
    DI = nan(nodeNum, nodeMax);
    coeff = nan(nodeNum, nodeMax);
    DIsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && nodeControl(i,j,1) == 0, continue; end
            if j>nodeNum && exControl(i,j-nodeNum,1) == 0, continue; end

            % autoregression plus other regression
            Yt = Y(1:(sigLen-p),i); % target
            Yti = ones(sigLen-p, p*2+1);
            for k=1:p
                Yti(:,1+k) = Y(k+1:sigLen-p+k,i); % target
                Yti(:,1+p+k) = Y(k+1:sigLen-p+k,j); % source
            end
            [b,bint,Yr] = regress(Yt,Yti);
            DIsub(i,j,1) = sum(b);
            DIsub(i,j,2) = sum(b(1:p+1));
            coeff(i,j) = DIsub(i,j,1)-DIsub(i,j,2); % actually this is sum of b(p+2:end)
            DI(i,j) = abs(coeff(i,j));
        end
    end
end

