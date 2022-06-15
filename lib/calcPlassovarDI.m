%%
% Calculate pLassoVAR (pairwise Lasso Vector Auto-Regression) DI
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
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end
    
    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc Pairwise LassoVAR
    DI = nan(nodeNum, nodeMax);
    coeff = nan(nodeNum, nodeMax);
    DIsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        [~,idx] = find(control(i,i,:)==1);
        Xt = Y(1:sigLen-lags,i);
        Yi = zeros(sigLen-lags, lags);
        for k=1:lags, Yi(:,k) = Y(1+k:sigLen-lags+k,i); end
        Xti = Yi(:,idx);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % autoregression plus other regression
            [~,idx] = find(control(i,j,:)==1);
            Yj = zeros(sigLen-lags, lags);
            for k=1:lags, Yj(:,k) = Y(1+k:sigLen-lags+k,j); end
            Xtj = [Xti, Yj(:,idx)];

            [b,~] = lasso(Xtj,Xt,'Lambda',lambda,'Alpha',elaAlpha);  % including Intercept
            DIsub(i,j,1) = sum(b);
            DIsub(i,j,2) = sum(b(1:size(Xti,2)));
            coeff(i,j) = DIsub(i,j,1)-DIsub(i,j,2); % actually this is sum of b(p+2:end)
            DI(i,j) = abs(coeff(i,j));
        end
    end
end

