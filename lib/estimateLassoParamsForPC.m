%%
% estimate Lasso (Elastic Net) Params for Partial Correlation
% returns lambda and alpha
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  subRate      subset sampling rate (0-1) (default:0.1)
%  cv           number of cross-validation (default:5)

function [lambda, elaAlpha, errMat] = estimateLassoParamsForPC(X, exSignal, nodeControl, exControl, subRate, cv, lambdaRange, alphaRange)
    if nargin < 6, cv = 5; end
    if nargin < 5, subRate = 0.1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    
    % get full training & test set for lasso regression
    trainFullSet = {};
    fullIdx = 1:nodeMax;
    for i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);

        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            x = Y(i,:).';
            y = Y(j,:).';
            idx = setdiff(nodeIdx,j);
            z = Y(idx,:).';
            trainFullSet{end+1} = {z, x};
            trainFullSet{end+1} = {z, y};
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
                    [x, z, testX, testZ] = getkFoldDataSetOfLassoPC(data{2}, data{1}, k, cv);

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

function [target, predict, testTarget, testPredict] = getkFoldDataSetOfLassoPC(orgTarget, orgPredict, k, N)
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
