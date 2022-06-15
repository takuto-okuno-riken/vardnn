%%
% Calculate Partial Correlation
% returns Partial Correlation (PC) and p-values (P)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC, P] = calcPartialCorrelation(X, exSignal, nodeControl, exControl, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    nodeMax = nodeNum + size(exSignal,1);
    
    % set node input
    Y = [X; exSignal];

    % check all same value or not
    for i=1:nodeMax
        A = unique(Y(i,:));
        if length(A)==1
            Y(i,mod(i,sigLen)+1) = A + 1.0e-8;
        end
    end

    if isempty(nodeControl) && isempty(exControl)
        [PC, P] = partialcorr(Y.');
    else
        PC = nan(nodeNum, nodeMax);
        P = nan(nodeNum, nodeMax);
        for i=1:nodeNum
            nodeIdx = [1:nodeNum];
            if ~isempty(nodeControl)
                [~,nodeIdx] = find(nodeControl(i,:)==1);
            end
            exIdx = [nodeNum+1:nodeNum+size(exSignal,1)];
            if ~isempty(exControl)
                [~,exIdx] = find(exControl(i,:)==1);
                exIdx = exIdx + nodeNum;
            end
            Y2 = Y([nodeIdx, exIdx],:);
            [PC2, P2] = partialcorr(Y2.');
            PC(i,[nodeIdx, exIdx]) = PC2(i,:);
            P(i,[nodeIdx, exIdx]) = P2(i,:);
        end
        for i=1:nodeNum
            for j=i:nodeNum, PC(j,i)=PC(i,j); P(j,i)=P(i,j); end
        end
    end

    % output control
    PC = PC(1:nodeNum,:); P = P(1:nodeNum,:); 
    if isFullNode == 0
        PC = PC(:,1:nodeNum); P = P(:,1:nodeNum); 
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        PC(:,1:nodeNum) = PC(:,1:nodeNum) .* nodeControl;
        P(:,1:nodeNum) = P(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        PC(:,nodeNum+1:end) = PC(:,nodeNum+1:end) .* exControl;
        P(:,nodeNum+1:end) = P(:,nodeNum+1:end) .* exControl;
    end
end
