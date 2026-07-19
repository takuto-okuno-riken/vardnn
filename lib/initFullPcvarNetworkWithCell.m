%%
% Calculate linear multivariate Vector Auto-Regression By PCA with Cell inputs of X and Ex.
% Only node connections have time lags, but not exogenous connections.
% This function assumes full connect control. exSignal only takes lag=1 (different from initMvarNetwork).
% returns mVAR network (net).
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  lags            number of lags for autoregression (default:3)
%  explainedTh     explained threshold for PCA components (default:0.99)
%  uniqueDecimal   taking unique value from conjuncted time series (option)
%  usehalf         use half for regress function (default:false)
%  usegpu          use gpu for regress function (default:false)
%  verbose         show verbose log (default:false)

function net = initFullPcvarNetworkWithCell(CX, CexSignal, lags, explainedTh, uniqueDecimal, usehalf, usegpu, verbose, usecache)
    if nargin < 9, usecache = false; end
    if nargin < 8, verbose = false; end
    if nargin < 7, usegpu = false; end
    if nargin < 6, usehalf = false; end
    if nargin < 5, uniqueDecimal = 0; end
    if nargin < 4, explainedTh = 0.99; end
    if nargin < 3, lags = 3; end
    if nargin < 2, CexSignal = {}; end
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    sigLen = size(CX{1},2);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    expTh = explainedTh * 100;

    % set vector auto-regression (VAR) inputs
    if verbose, disp('set VAR inputs'); end
    allInLen = 0;
    for i=1:cxNum
        allInLen = allInLen + size(CX{i},2) - lags;
    end
    Y = single(nan(allInLen,nodeNum));
    Xti = single(nan(allInLen,nodeNum*lags));
    Xe = single(nan(allInLen,exNum));
    xts = 1;
    for i=1:cxNum
        % set node input
        X = flipud(CX{i}.'); % need to flip signal

        sLen = size(X,1);
        sl = sLen-lags;
        Xt = single(zeros(sl, lags*nodeNum));
        for p=1:lags
            Xt(:,1+nodeNum*(p-1):nodeNum*p) = X(1+p:sl+p,:);
        end
        Y(xts:xts+sl-1,:) = X(1:sl,:);
        Xti(xts:xts+sl-1,:) = Xt;
        if exNum > 0
            E = flipud(CexSignal{i}.');
            Xe(xts:xts+sl-1,:) = E(1:sl,:);
        end
        xts = xts + sl;
    end
    if uniqueDecimal > 0
        Y = single(int32(Y / uniqueDecimal)) * uniqueDecimal;
        Xti = single(int32(Xti / uniqueDecimal)) * uniqueDecimal;
        Xe = single(int32(Xe / uniqueDecimal)) * uniqueDecimal;
    end

    clear X; % memory clear
    clear Xt; % memory clear

    % calculate mean and covariance of each node
    Z = [];
    for i=1:cxNum, Z = [Z, CX{i}]; end
    if usegpu
        % check gpu maxGridSize
        [sz1,sz2] = size(Z);
        maxGridSize = gpuDevice().MaxGridSize(1);
        if sz1*sz2 > maxGridSize
            disp(['error : input matrix size exceeded gpu device MaxGridSize : ' num2str(maxGridSize)]);
            disp('gpu device will not be used.');
            usegpu = false;
        else
            Z = gpuArray(single(Z));
        end
    end
    cxM = mean(Z.');
    cxCov = cov(Z.',1);
    clear Z; % memory clear

    clear CX; % memory clear
    clear CexSignal; % memory clear

    % apply the Principal Component Regress function
    if verbose, disp('apply PCA'); end
    cacheName = ['results/glm/initglmpc-pca-cache-' num2str(size(Xti,1)) 'x' num2str(size(Xti,2)) '-l' num2str(lags) '-d' num2str(uniqueDecimal) '-h' num2str(usehalf) '.mat'];
    if usegpu
        availableMemory = gpuDevice().AvailableMemory / 4;
        if availableMemory > size(Xti,1) * size(Xti,2) * 8
            Xti = gpuArray(Xti);
        end
    end
    if ~usegpu && usecache && exist(cacheName,'file') 
        load(cacheName);
    else
        [~,score,~,~,explained] = pca(double(Xti)); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        if ~usegpu && usecache
            save(cacheName,'score','explained','-v7.3');
        end
    end

    % find 99% component range
    expTotal = 0;
    maxComp = size(score,2);
    for j=1:size(Xti,2)
        expTotal = expTotal + explained(j);
        if expTotal >= expTh
            maxComp = j;
            break;
        end
    end

    % preparation
    if verbose, disp('prepare regression inputs'); end
    score = score(:,1:maxComp); % reduce memory
    if exNum > 0
        pcXti = [score, Xe, ones(allInLen,1)];
        X = [Xti, Xe, ones(allInLen,1)];
    else
        pcXti = [score, ones(allInLen,1)];
        X = [Xti, ones(allInLen,1)];
    end
    Xtiisc = invQR(Xti) * score;
    clear score;
    clear explained;
    clear Xe; % memory clear
    clear Xti;

    if usegpu
        pcXti = gpuArray(pcXti);
        Y = gpuArray(Y);
        X = gpuArray(X);
    end
    [~, ~, perm, RiQ, dR2i] = regressPrepare(pcXti);
    df = size(X,1) - size(X,2); % degree of freedom
    if df < 1
        df = 1;
        dX2i = diag(inv(double(X'*X))); % this will be pseudo inv
    else
        dX2i = diag(invQR(X'*X)); % get inv(X'*X)
    end
    clear X;

    if usehalf && ~usegpu
        pcXti = half(pcXti);
        Y = half(Y);
    end

    % apply the regress function
    b = cell(nodeNum,1);
    r = cell(nodeNum,1);
    T = cell(nodeNum,1);

%    for n=1:nodeNum
%    if isempty(gcp('nocreate'))
%        parpool('Threads');   % this doesn't work well. perhaps, threads are used in chol, inv, etc functions. multi-process is better.
%    end
    parfor n=1:nodeNum
        if verbose, disp(['calc node' num2str(n)]); end

        Xt = Y(:,n);
        [a, r{n}] = regressLinear(Xt, pcXti, [], [], perm, RiQ, dR2i);
        bt = Xtiisc * a(1:maxComp); % last one is coefficient of intercept
        pr = single(r{n});
        s  = sqrt(sum(pr.*pr)/df); pr = []; % mem clear
        se = sqrt(dX2i * (s*s));
        b2 = [];
        s2 = [];
        for p=1:lags
            if p==1
                b2 = [b2; bt(1+(p-1)*nodeNum:p*nodeNum); a(1+maxComp:end-1)];
                s2 = [s2; se(1+(p-1)*nodeNum:p*nodeNum); se(1+lags*nodeNum:end-1)];
            else
                b2 = [b2; bt(1+(p-1)*nodeNum:p*nodeNum); zeros(exNum,1)];
                s2 = [s2; se(1+(p-1)*nodeNum:p*nodeNum); ones(exNum,1)];
            end
        end
        s2 = [s2; se(end)];
        b{n} = [b2; a(end)]; % last one is coefficient of intercept
        T{n}  = b{n}./s2;
        if usegpu
            b{n} = gather(b{n});
            r{n} = gather(r{n});
            T{n} = gather(T{n});
        end
        if usehalf 
            b{n} = half(b{n});
            r{n} = half(r{n});
            T{n} = half(T{n});
        end
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.sigLen = sigLen;
    net.cxM = cxM;
    net.cxCov = cxCov;
    net.lags = lags;
    net.bvec = b;
    net.rvec = r;
    net.Tvec = T;
end
