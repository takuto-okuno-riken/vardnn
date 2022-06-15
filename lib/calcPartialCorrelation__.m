%%
% Calculate Partial Correlation by linear regression
% returns Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = calcPartialCorrelation__(X, exSignal, nodeControl, exControl, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    
    PC = nan(nodeNum,nodeMax);
    fullIdx = 1:nodeMax;
    parfor i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);
        x = Y(i,:).';
        
        A = nan(1,nodeMax,class(X));
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            idx = setdiff(nodeIdx,j);
            z = [Y(idx,:).', ones(sigLen,1)]; % add intercept

            [~, ~, perm, RiQ] = regressPrepare(z);
            [~, r1] = regressLinear(x, z, [], [], perm, RiQ);
            [~, r2] = regressLinear(Y(j,:).', z, [], [], perm, RiQ);

            A(j) = (r1.'*r2) / (sqrt(r1.'*r1)*sqrt(r2.'*r2));
        end
        PC(i,:) = A;
    end
    A = ones(nodeNum,'logical'); A = tril(A,-1);
    idx = find(A==1);
    B = PC(:,1:nodeNum); C = B';
    B(idx) = C(idx);
    PC(:,1:nodeNum) = B;

    % output control
    PC = PC(1:nodeNum,:);
    if isFullNode == 0
        PC = PC(:,1:nodeNum);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        PC(:,1:nodeNum) = PC(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        PC(:,nodeNum+1:end) = PC(:,nodeNum+1:end) .* exControl;
    end
end
