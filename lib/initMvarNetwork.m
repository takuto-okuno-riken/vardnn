%%
% calculate linear multivariate vector auto-regression weights and Create MVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)

function net = initMvarNetwork(X, exSignal, nodeControl, exControl, lags)
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
    bint = cell(nodeNum,1);
    r = cell(nodeNum,1);
    rint = cell(nodeNum,1);
    stats = cell(nodeNum,1);

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
        exIdx = [nodeNum+1:nodeNum+size(exSignal,1)];
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
        Xti = [Yj(:,idx), ones(sigLen-p,1)]; % might not be good to add bias
        % apply the regress function
        [b{i},bint{i},r{i},rint{i},stats{i}] = regress(Xt,Xti);
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.bvec = b;
    net.bint = bint;
    net.rvec = r;
    net.rint = rint;
    net.stats = stats;
end