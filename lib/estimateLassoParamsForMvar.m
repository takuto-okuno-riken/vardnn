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
    inputNum = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end

    % get full training & test set for lasso regression
    trainFullSet = {};
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);
        
        % PLS vector auto-regression (VAR)
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);
        trainFullSet{end+1} = {Xti, Xt};

        for j=1:inputNum
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % lasso vector auto-regression (VAR)
            control2 = control;
            control2(i,j,:) = 0;
            [~,idx2] = find(control2(i,:,:)==1);
            Xtj = Yj(:,idx2);

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
            for n=1:length(trainSet)
                for k=1:cv
                    data = trainSet{n};
                    [x, z, testX, testZ] = getkFoldDataSetOfLassoGC(data{2}, data{1}, k, cv);

                    [b,info] = lasso(z,x,'Lambda',lambda,'Alpha',elaAlpha); % including Intercept
                    r = x - (z*b + info.Intercept);
                    r1 = [r1; r];
                    r = testX - (testZ*b + info.Intercept);
                    r2 = [r2; r];
                end
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
