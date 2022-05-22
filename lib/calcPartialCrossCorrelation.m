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
        Ulen(i) = length(unique(Y(i,:)));
    end
    
    NPCC = nan(nodeNum,nodeMax,maxlag*2+1);
    fullIdx = 1:nodeMax;
%    for i=1:nodeNum
    parfor i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);
        
        A = nan(nodeMax,maxlag*2+1);
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

                [Q, R, perm, RiQ] = regressPrepare(z);
                [~, r1] = regressLinear(x, z, Q, R, perm, RiQ);
                [~, r2] = regressLinear(y, z, Q, R, perm, RiQ);

                [A(j,:), ~] = xcov(r1,r2,maxlag,'normalized');
            end
        end
        NPCC(i,:,:) = A;
    end
    for i=1:nodeNum
        for j=i:nodeMax, NPCC(j,i,:) = NPCC(i,j,:); end
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
