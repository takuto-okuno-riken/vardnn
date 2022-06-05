%%
% Caluclate Normalized Partial Cross-Correlation
% returns Normalized Partial Cross-Correlation (NPCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [NPCC, lags] = calcPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, isFullNode)
    if nargin < 6, isFullNode = 0; end
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
    for i=1:nodeMax
        Ulen(i) = length(unique(single(Y(i,:)))); % 'half' does not support
    end
    
    NPCC = nan(nodeNum,nodeMax,maxlag*2+1,class(X));
    fullIdx = 1:nodeMax;
%    for i=1:nodeNum
    parfor i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);
        
        A = nan(nodeMax,maxlag*2+1,class(X));
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            if Ulen(i)==1 && Ulen(j)==1
                A(j,:) = 0;
            else
                x = Y(i,:).';
                y = Y(j,:).';
                idx = setdiff(nodeIdx,j);
                z = [Y(idx,:).', ones(sigLen,1)];

                [~, ~, perm, RiQ] = regressPrepare(z);
                [~, r1] = regressLinear(x, z, [], [], perm, RiQ);
                [~, r2] = regressLinear(y, z, [], [], perm, RiQ);

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
