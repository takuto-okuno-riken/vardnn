%%
% calculate linear pairwise Principal Component Vector Auto-Regression weights and Create pPCVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  explainedTh     explained threshold for PCA components (default:0.99)

function net = initPpcvarNetwork(X, exSignal, nodeControl, exControl, lags, explainedTh)
    if nargin < 6, explainedTh = 0.99; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    expTh = explainedTh * 100;
    p = lags;

    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;
    
    coeff = cell(nodeNum,nodeMax);
    score = cell(nodeNum,nodeMax);
    latent = cell(nodeNum,nodeMax);
    explained = cell(nodeNum,nodeMax);
    mu = cell(nodeNum,nodeMax);
    maxComp = cell(nodeNum,nodeMax);
    b = cell(nodeNum,nodeMax);
    bint = cell(nodeNum,nodeMax);
    r = cell(nodeNum,nodeMax);
    rint = cell(nodeNum,nodeMax);
    stats = cell(nodeNum,nodeMax);

    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            Y1 = flipud(Y(i,:));
            Y2 = flipud(Y(j,:));

            % autoregression plus other regression
            Yt = Y2(1:sigLen-p).'; % TODO: X1 & X2 opposite ??
            Yti = ones(sigLen-p, p*2);
            for k=1:p
                Yti(:,k) = Y2(k+1:sigLen-p+k);
                Yti(:,p+k) = Y1(k+1:sigLen-p+k);
            end

            % apply the Principal Component Regress function
            [coeff{i,j},score{i,j},latent{i,j},~,explained{i,j},mu{i,j}] = pca(Yti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

            % find 99% component range
            expTotal = 0;
            for k=1:size(Yti,2)
                expTotal = expTotal + explained{i,j}(k);
                if expTotal > expTh
                    maxComp{i,j} = k;
                    break;
                end
            end
            pcXti = [score{i,j}(:,1:maxComp{i,j}), ones(sigLen-p,1)]; % might not be good to add bias
            [b{i,j},bint{i,j},r{i,j},rint{i,j},stats{i,j}] = regress(Yt, pcXti);
        end
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.coeff = coeff;
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