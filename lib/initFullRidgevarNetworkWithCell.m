%%
% calculate linear multivariate Ridge Vector Auto-Regression with Cell inputs of X and Ex.
% Only node connections have time lags, but not exogenous connections.
% This function assumes full connect control.
% returns mVAR network (net). exSignal only takes lag=1 (different from initMvarNetwork).
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  lags            number of lags for autoregression (default:3)
%  k               k for ridge regression (default:0.1)
%  uniqueDecimal   taking unique value from conjuncted time series (option)
%  usehalf         use half for regress function (default:false)
%  usegpu          use gpu for regress function (default:false)
%  verbose         show verbose log (default:false)

function net = initFullRidgevarNetworkWithCell(CX, CexSignal, lags, k, uniqueDecimal, usehalf, usegpu, verbose)
    if nargin < 8, verbose = false; end
    if nargin < 7, usegpu = false; end
    if nargin < 6, usehalf = false; end
    if nargin < 5, uniqueDecimal = 0; end
    if nargin < 4, k = 0.1; end
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
        Xti = [Xti, Xe, ones(allInLen,1)];
    else
        Xti = [Xti, ones(allInLen,1)];
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
        % matlab ridge() compatible version
        Cp = [Xti; k * eye(size(Xti,2))];
        Cp2 = double(Cp' * Cp);
        dc = det(Cp2);
        if dc == 0 || isinf(dc)
            Cp2i = invQR(Cp2);
        else
            Cp2i = inv(Cp2);
        end
        Cp2iCp = Cp2i * Cp';
        clear Cp; clear Cp2; clear Cp2i;

        % for T-value based on Cule et al., 2011
        A = double(Xti' * Xti);
        B = A + k * eye(size(Xti,2));
        dc = det(B);
        if dc == 0 || isinf(dc)
            Bi = invQR(B);
        else
            Bi = inv(B);
        end
        dX2i = diag(Bi*A*Bi);
        clear A; clear B; clear Bi;

        if ~usegpu && usecache
            save(cacheName,'Cp2iCp','-v7.3');
        end
    end
    if usehalf && ~usegpu
        Xti = half(Xti);
        Y = half(Y);
    end

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
        bt =  single(Cp2iCp * [y; zeros(size(Xti,2),1)]); % matlab ridge() compatible version
%        bt2 = [Xti; k * eye(size(Xti,2))] \ [y; zeros(size(Xti,2),1)]; % ridge() compatible version. this is slow.
        r{n} = y - Xti*bt;
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
        b{n} = [b2; bt(end)];
        s2 = [s2; se(end)];
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
