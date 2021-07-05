%%
% estimate Lasso (Elastic Net) Params for mVAR
% returns lambda and alpha
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  subRate      subset sampling rate (0-1) (default:0.1)
%  cv           number of cross-validation (default:5)

function [lambda, elaAlpha, errMat] = estimateLassoParamsForMvar(X, exSignal, nodeControl, exControl, lags, subRate, cv, lambdaRange, alphaRange)
    if nargin < 6, cv = 5; end
    if nargin < 5, subRate = 0.1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;
    p = lags;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-p, p*nodeMax);
    for k=1:p
        Yj(:,1+nodeMax*(k-1):nodeMax*k) = Y(1+k:sigLen-p+k,:);
    end

    % get full training & test set for lasso regression
    trainFullSet = {};
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeMax];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:p
            idx = [idx, nodeIdx+nodeMax*(k-1), exIdx+nodeMax*(k-1)];
        end
        idxList = [nodeIdx, exIdx];
        nlen = length(idxList);
        Xt = Y(1:sigLen-p,i);
        Xti = Yj(:,idx);
        trainFullSet{end+1} = {Xti, Xt};

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            % PLS vector auto-regression (VAR)
            jIdx = idx;
            for k=p:-1:1
                jIdx(j+nlen*(k-1)) = [];
            end
            Xtj = Yj(:,jIdx);
            trainFullSet{end+1} = {Xtj, Xt};
        end
    end
    
    % get subset of training & test set for lasso regression
    tlen = length(trainFullSet);
    idx = randperm(tlen,floor(tlen*subRate));
    trainSet = trainFullSet(idx);

    % cross validation
    llen = length(lambdaRange);
    alen = length(alphaRange);
    errMat = nan(llen,alen);
    r1All = cell(llen,alen);
    r2All = cell(llen,alen);
    minErr = 100;
    minPair = [];
    for i=1:llen
        for j=1:alen
            lambda = lambdaRange(i);
            elaAlpha = alphaRange(j);
            
            % cross validation
            r1 = [];
            r2 = [];
            for k=1:cv
                data = trainSet{k};
                [x, z, testX, testZ] = getkFoldDataSetOfLassoGC(data{2}, data{1}, k, cv);

                [b,info] = lasso(z,x,'Lambda',lambda,'Alpha',elaAlpha); % including Intercept
                r = x - (z*b + info.Intercept);
                r1 = [r1; r];
                r = testX - (testZ*b + info.Intercept);
                r2 = [r2; r];
            end
            rmse = sqrt(r2.'*r2 / length(r2));
            errMat(i,j) = rmse;
            r1All{i,j} = r1;
            r2All{i,j} = r2;
            if minErr > rmse
                minErr = rmse;
                minPair = [i,j];
            end
        end
    end
    lambda = lambdaRange(minPair(1));
    elaAlpha = alphaRange(minPair(2));
end

function [target, predict, testTarget, testPredict] = getkFoldDataSetOfLassoGC(orgTarget, orgPredict, k, N)
    sampleNum = size(orgTarget,1);
    un = floor(sampleNum / N);
    st = (k-1)*un+1;
    ed = k*un;
    if k==N, ed = sampleNum; end

    fullIdx = 1:sampleNum;
    idx = setdiff(fullIdx,st:ed);
    target = orgTarget(idx,:);
    predict = orgPredict(idx,:);
    testTarget = orgTarget(st:ed,:);
    testPredict = orgPredict(st:ed,:);
end
