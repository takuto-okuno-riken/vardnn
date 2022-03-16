%%
% calculate linear multivariate Principal Component Vector Auto-Regression weights and Create mPCVAR network
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  explainedTh     explained threshold for PCA components (default:0.99)
%  uniqueDecimal   taking unique value from conjuncted time series (option)

function net = initMpcvarNetworkWithCell(CX, CexSignal, nodeControl, exControl, lags, explainedTh, uniqueDecimal)
    if nargin < 7, uniqueDecimal = 0; end
    if nargin < 6, explainedTh = 0.99; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, CexSignal = {}; end
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    sigLen = size(CX{1},2);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    inputNum = nodeNum + exNum;
    expTh = explainedTh * 100;

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    coeff = cell(nodeNum,1);
%    score = cell(nodeNum,1);
    latent = cell(nodeNum,1);
    explained = cell(nodeNum,1);
    mu = cell(nodeNum,1);
    maxComp = cell(nodeNum,1);
    b = cell(nodeNum,1);
    r = cell(nodeNum,1);
    stats = cell(nodeNum,1);
    
    % set vector auto-regression (VAR) inputs
    idxs = cell(nodeNum,1);
    for n=1:nodeNum
        [~,idxs{n}] = find(control(n,:,:)==1);
    end
    allInLen = 0;
    for i=1:cxNum
        allInLen = allInLen + size(CX{i},2) - lags;
    end

    % calculate mean and covariance of each node
    Y = [];
    for i=1:cxNum, Y = [Y, CX{i}]; end
    cxM = mean(Y.');
    cxCov = cov(Y.');
    Y = []; % memory clear

    % apply the Principal Component Regress function
%    for n=1:nodeNum
    parfor n=1:nodeNum
        Xt = single(nan(allInLen,1));
        Xti = single(nan(allInLen,length(idxs{n})));
        xts = 1;

        % this is redundant for every node. but it is necessary to avoid
        % too much memory consumption
        for i=1:cxNum
            % set node input
            Y = single(CX{i});
            if exNum > 0
                Y = [Y; single(CexSignal{i})];
            end
            Y = flipud(Y.'); % need to flip signal

            sLen = size(Y,1);
            sl = sLen-lags;
            Yj = single(zeros(sl, lags*inputNum));
            for p=1:lags
                Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sl+p,:);
            end
            Xt(xts:xts+sl-1,:) = Y(1:sl,n);
            Xti(xts:xts+sl-1,:) = Yj(:,idxs{n});
            xts = xts + sl;
        end
        Y = [];  % clear memory
        Yj = []; % clear memory
        
        if uniqueDecimal > 0
            A = int32([Xt, Xti] / uniqueDecimal);
            A = unique(A,'rows');
            Xt = single(A(:,1)) * uniqueDecimal;
            Xti = single(A(:,2:end)) * uniqueDecimal;
            A = []; % clear memory
        end

        [coeff{n},score,latent{n},~,explained{n},mu{n}] = pca(Xti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

        % find 99% component range
        expTotal = 0;
        maxComp{n} = size(score,2);
        for j=1:size(Xti,2)
            expTotal = expTotal + explained{n}(j);
            if expTotal >= expTh
                maxComp{n} = j;
                break;
            end
        end
        pcXti = [score(:,1:maxComp{n}), ones(size(Xti,1),1)]; % might not be good to add bias
        [b{n},~,r{n},~,stats{n}] = regress(Xt, pcXti);
        b{n} = single(b{n});
        pcXti = [];  % clear memory
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.sigLen = sigLen;
    net.cxM = cxM;
    net.cxCov = cxCov;
    net.lags = lags;
    net.coeff = coeff;
    net.latent = latent;
    net.explained = explained;
    net.explainedTh = explainedTh;
    net.mu = mu;
    net.maxComp = maxComp;
    net.bvec = b;
    net.rvec = r;
    net.stats = stats;
end
