%%
% Caluclate Partial Correlation
% returns Partial Correlation (PC) and p-values (P)
% input:
%  X     multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC, P] = calcPartialCorrelation(X, exSignal, nodeControl, exControl, isFullNode)
    if nargin < 5
        isFullNode = 0;
    end
    if nargin < 4
        exControl = [];
    end
    if nargin < 3
        nodeControl = [];
    end
    if nargin < 2
        exSignal = [];
    end
    nodeNum = size(X,1);
    
    % set node input
    if ~isempty(exSignal) && isFullNode > 0
        X = [X; exSignal];
    end
    [PC, P] = partialcorr(X.');

    if ~isempty(exSignal)
        PC = PC(1:nodeNum,:); P = P(1:nodeNum,:); 
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
