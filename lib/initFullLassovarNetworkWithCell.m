%%
% calculate linear multivariate Lasso Vector Auto-Regression with Cell inputs of X and Ex.
% Only node connections have time lags, but not exogenous connections.
% This function assumes full connect control. exSignal only takes lag=1 (different from initMvarNetwork).
% returns mVAR network (net).
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  lags            number of lags for autoregression (default:3)
%  DFrate          nonzero B value rate for lasso regression (default:0.1)
%  uniqueDecimal   taking unique value from conjuncted time series (option)
%  usehalf         use half for regress function (default:false)
%  usegpu          use gpu for regress function (default:false)
%  verbose         show verbose log (default:false)

function net = initFullLassovarNetworkWithCell(CX, CexSignal, lags, DFrate, uniqueDecimal, usehalf, usegpu, verbose)
    if nargin < 8, verbose = false; end
    if nargin < 7, usegpu = false; end
    if nargin < 6, usehalf = false; end
    if nargin < 5, uniqueDecimal = 0; end
    if nargin < 4, DFrate = 0.1; end
    if nargin < 3, lags = 3; end
    if nargin < 2, CexSignal = {}; end
    usecache = false;
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    sigLen = size(CX{1},2);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end

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
        y = single(zeros(sl, lags*nodeNum));
        for p=1:lags
            y(:,1+nodeNum*(p-1):nodeNum*p) = X(1+p:sl+p,:);
        end
        Y(xts:xts+sl-1,:) = X(1:sl,:);
        Xti(xts:xts+sl-1,:) = y;
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
    clear y; % memory clear

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

    % preparation
    if verbose, disp('prepare regression inputs'); end
    if exNum > 0
        Xti = [Xti, Xe];
    else
        Xti = [Xti];
    end
    clear Xe; % memory clear

    if usegpu
        Xti = gpuArray(Xti);
        Y = gpuArray(Y);
    end
    cacheName = ['results/glm/initglm-rdgregpre-cache-' num2str(size(Xti,1)) 'x' num2str(size(Xti,2)) '-l' num2str(lags) '-d' num2str(uniqueDecimal) '-h' num2str(usehalf) '.mat'];
    if ~usegpu && usecache && exist(cacheName,'file') 
        load(cacheName);
    else
        % for T-value
        A = double(Xti' * Xti);
        dc = det(A);
        if dc == 0 || isinf(dc)
            Ai = invQR(A);
        else
            Ai = inv(A);
        end
        dX2i = diag(Ai);
        clear A; clear Ai;

        if ~usegpu && usecache
            save(cacheName,'dX2i','-v7.3');
        end
    end
    if usehalf && ~usegpu
        Xti = half(Xti);
        Y = half(Y);
    end

    % DF num for DFmax
    DFnum = DFrate * size(Xti,2);

    % apply the regress function
    b = cell(nodeNum,1);
    r = cell(nodeNum,1);
    T = cell(nodeNum,1);
    df = size(Xti,1) - size(Xti,2); % degree of freedom

%    for n=1:nodeNum
%    if isempty(gcp('nocreate'))
%        parpool('Threads');   % this doesn't work well. perhaps, threads are used in chol, inv, etc functions. multi-process is better.
%    end
    parfor n=1:nodeNum
        if verbose, disp(['calc node' num2str(n)]); end

        y = Y(:,n);
        [bt, stat] = lasso(Xti, y, 'DFmax', DFnum, 'NumLambda', 10); % sampling lambda test
        DFs = stat.DF - DFnum;
        [d,idx] = min(abs(DFs));
        bt = bt(:,idx);
        it = stat.Intercept(idx);
        r{n} = y - Xti*bt - it;
        pr = single(r{n});
        s  = sqrt(sum(pr.*pr)/df);
        se = sqrt(dX2i * (s*s));
        % reconstruct for compatibility
        b2 = [];
        s2 = [];
        for p=1:lags
            if exNum > 0
                if p==1
                    b2 = [b2; bt(1+nodeNum*(p-1):nodeNum*p); bt(1+nodeNum*lags:end-1)];
                    s2 = [s2; se(1+(p-1)*nodeNum:p*nodeNum); se(1+lags*nodeNum:end-1)];
                else
                    b2 = [b2; bt(1+nodeNum*(p-1):nodeNum*p); zeros(exNum,1)];
                    s2 = [s2; se(1+(p-1)*nodeNum:p*nodeNum); ones(exNum,1)];
                end
            else
                b2 = [b2; bt(1+(p-1)*nodeNum:p*nodeNum)];
                s2 = [s2; se(1+(p-1)*nodeNum:p*nodeNum)];
            end
        end
        b{n} = [b2; it];
        s2 = [s2; Inf]; % this case, intercept has dummy T-value
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
