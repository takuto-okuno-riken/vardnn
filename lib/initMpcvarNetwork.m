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
    inputNum = nodeNum + exNum;
    
    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    Y = flipud(Y.'); % need to flip signal
    
    coeff = cell(nodeNum,1);
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

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);

        % vector auto-regression (VAR)
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);

        % apply the Principal Component Regress function
        [coeff{i},score{i},latent{i},~,explained{i},mu{i}] = pca(Xti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

        % find 99% component range
        expTotal = 0;
        maxComp{i} = size(score,2);
        for j=1:size(Xti,2)
            expTotal = expTotal + explained{i}(j);
            if expTotal >= expTh
                maxComp{i} = j;
                break;
            end
        end
        pcXti = [score{i}(:,1:maxComp{i}), ones(sigLen-lags,1)]; % might not be good to add bias
        [b{i},bint{i},r{i},rint{i},stats{i}] = regress(Xt, pcXti);
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.coeff = coeff;
    net.score = score;
    net.latent = latent;
    net.explained = explained;
    net.explainedTh = explainedTh;
    net.mu = mu;
    net.maxComp = maxComp;
    net.bvec = b;
    net.bint = bint;
    net.rvec = r;
    net.rint = rint;
    net.stats = stats;
end