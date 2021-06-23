%%
% calculate linear multivariate Lasso Vector Auto-Regression weights and Create mLassoVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  lambda          lambda for the Lasso (default:0.01)
%  elaAlpha        Elastic Net Alpha for the Lasso (default:1)

function net = initMlassovarNetwork(X, exSignal, nodeControl, exControl, lags, lambda, elaAlpha)
    if nargin < 7, elaAlpha = 1; end
    if nargin < 6, lambda = 0.01; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;
    
    b = cell(nodeNum,1);
    info = cell(nodeNum,1);

    p = lags;
    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-p, p*nodeMax);
    for k=1:p
        Yj(:,1+nodeMax*(k-1):nodeMax*k) = Y(1+k:sigLen-p+k,:);
    end
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeNum+exNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:p
            idx = [idx, nodeIdx+nodeMax*(k-1), exIdx+nodeMax*(k-1)];
        end

        % vector auto-regression (VAR)
        Xt = Y(1:sigLen-p,i);
        Xti = Yj(:,idx);
        % apply the regress function
        [b{i}, info{i}] = lasso(Xti,Xt,'Lambda',lambda,'Alpha',elaAlpha); % including Intercept
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.lambda = lambda;
    net.elaAlpha = elaAlpha;
    net.bvec = b;
    net.info = info;
end