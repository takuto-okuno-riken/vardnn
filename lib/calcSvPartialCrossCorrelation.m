%%
% Caluclate Support Vector Normalized Partial Cross-Correlation (SvNPCC)
% returns Support Vector Normalized Partial Cross-Correlation (NPCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  kernel       kernel for SVM (default:'linear', 'gaussian', 'rbf')
%  kernelScale  kernelScale for SVM (default:'auto', 1)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [NPCC, lags] = calcSvPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, kernel, kernelScale, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, kernelScale = 'auto'; end
    if nargin < 6, kernel = 'linear'; end
    if nargin < 5, maxlag = 5; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];

    % check all same value or not
    uLen = zeros(nodeMax,1);
    for i=1:nodeMax
        uLen(i) = length(unique(single(Y(i,:)))); % 'half' does not support
    end
    
    % each pool have each own gloval value.
%   setSvPCCFlag(false); % init global value for 'for loop'
%%{
    p = gcp; % If no pool, create new
    parfor i=1:p.NumWorkers
        setSvPCCFlag(false); % init global value for 'parfor'
    end
%%}
    NPCC = nan(nodeNum,nodeMax,maxlag*2+1,class(X));
    fullIdx = 1:nodeMax;
%    for i=1:nodeNum
    parfor i=1:nodeNum
        if getSvPCCFlag() == true, continue; end
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);

        A = nan(nodeMax,maxlag*2+1,class(X));
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            if getSvPCCFlag() == true, continue; end

            if uLen(i)==1 && uLen(j)==1
                A(j,:) = 0;
            else
                x = single(Y(i,:).'); % 'half' does not support
                y = single(Y(j,:).'); % 'half' does not support
                idx = setdiff(nodeIdx,j);
                z = single(Y(idx,:).'); % 'half' does not support

                mdl = fitrsvm(z,x,'KernelFunction',kernel,'KernelScale',kernelScale); %,'Standardize',true); % bias will be calcurated
                Si = predict(mdl, z);
                r1 = Si - x;

                % if r1 are all NaN, finish this function
                if length(x) == sum(isnan(r1))
                    setSvPCCFlag(true); continue;
                end

                mdl = fitrsvm(z,y,'KernelFunction',kernel,'KernelScale',kernelScale); %,'Standardize',true); % bias will be calcurated
                Si = predict(mdl, z);
                r2 = Si - y;

                [A(j,:), ~] = xcov(r1,r2,maxlag,'normalized');
            end
        end
        NPCC(i,:,:) = A;
    end
    A = ones(nodeNum,'logical'); A = tril(A,-1);
    idx = find(A==1);
    for i=1:size(NPCC,3)
        B = NPCC(:,1:nodeNum,i); C = B';
        B(idx) = C(idx);
        NPCC(:,1:nodeNum,i) = B;
    end
    lags = -maxlag:maxlag;

    % output control
    NPCC = NPCC(1:nodeNum,:,:);
    if isFullNode == 0
        NPCC = NPCC(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        NPCC(:,1:nodeNum,:) = NPCC(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        NPCC(:,nodeNum+1:end,:) = NPCC(:,nodeNum+1:end,:) .* exControl;
    end
end

function setSvPCCFlag(val)
    global svPCCFlag;
    svPCCFlag = val;
%    disp(['set svPCFlag = ' num2str(val)]);
end

function val = getSvPCCFlag()
    global svPCCFlag;
    val = svPCCFlag;
%    disp(['get svPCFlag = ' num2str(val)]);
end
