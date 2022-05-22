%%
% Caluclate Random Forest Partial Correlation
% returns Random Forest Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  numTrees     number of trees for random forest (default:10)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = calcRfPartialCorrelation(X, exSignal, nodeControl, exControl, numTrees, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, numTrees = 10; end
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
        A = unique(Y(i,:));
        if length(A)==1
            Y(i,mod(i,sigLen)+1) = A + 1.0e-8;
        end
    end

    fullIdx = 1:nodeMax;
    PC = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);

        for j=i:nodeMax
%        parfor j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            x = Y(i,:).';
            y = Y(j,:).';
            idx = setdiff(nodeIdx,j);
            z = Y(idx,:).';

            mdl = TreeBagger(numTrees,z,x,'Method','regression');
            Si = predict(mdl, z);
            r1 = Si - x;

            mdl = TreeBagger(numTrees,z,y,'Method','regression');
            Si = predict(mdl, z);
            r2 = Si - y;

            PC(i,j) = (r1.'*r2) / (sqrt(r1.'*r1)*sqrt(r2.'*r2));
        end
    end
    for i=1:nodeNum
        for j=i:nodeMax, PC(j,i) = PC(i,j); end
    end

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
