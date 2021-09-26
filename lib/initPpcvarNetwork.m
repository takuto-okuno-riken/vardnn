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
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

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

            % apply the Principal Component Regress function
            [coeff{i,j},score{i,j},latent{i,j},~,explained{i,j},mu{i,j}] = pca(Xtj); % relation : Xtj == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

            % find 99% component range
            expTotal = 0;
            maxComp{i,j} = size(Xtj,2);
            for k=1:size(Xtj,2)
                expTotal = expTotal + explained{i,j}(k);
                if expTotal >= expTh
                    maxComp{i,j} = k;
                    break;
                end
            end
            pcXtj = [score{i,j}(:,1:maxComp{i,j}), ones(sigLen-lags,1)]; % might not be good to add bias
            [b{i,j},bint{i,j},r{i,j},rint{i,j},stats{i,j}] = regress(Xt, pcXtj);
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