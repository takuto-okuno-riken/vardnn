%%
% calculate linear multivariate Principal Component Vector Auto-Regression weights and Create mPCVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  explainedTh     explained threshold for PCA components (default:0.99)

function net = initMpcvarNetwork(X, exSignal, nodeControl, exControl, lags, explainedTh)
    if nargin < 6, explainedTh = 0.99; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    expTh = explainedTh * 100;

    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;
    
    coeff = cell(nodeNum,1);
    invcoeff = cell(nodeNum,1);
    score = cell(nodeNum,1);
    latent = cell(nodeNum,1);
    explained = cell(nodeNum,1);
    mu = cell(nodeNum,1);
    maxComp = cell(nodeNum,1);
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
        Xti = Yj(:,idx);

        % apply the Principal Component Regress function
        [coeff{i},score{i},latent{i},~,explained{i},mu{i}] = pca(Xti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        invcoeff{i} = inv(coeff{i}); % this does full size inversion of matrix. might not be good.

        % find 99% component range
        expTotal = 0;
        for j=1:size(Xti,2)
            expTotal = expTotal + explained{i}(j);
            if expTotal > expTh
                maxComp{i} = j;
                break;
            end
        end
        pcXti = [score{i}(:,1:maxComp{i}), ones(sigLen-p,1)]; % might not be good to add bias
        [b{i},bint{i},r{i},rint{i},stats{i}] = regress(Xt, pcXti);
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.coeff = coeff;
    net.invcoeff = invcoeff;
    net.score = score;
    net.latent = latent;
    net.explained = explained;
    net.mu = mu;
    net.maxComp = maxComp;
    net.bvec = b;
    net.bint = bint;
    net.rvec = r;
    net.rint = rint;
    net.stats = stats;
end