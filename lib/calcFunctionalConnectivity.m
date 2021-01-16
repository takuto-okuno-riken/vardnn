%%
% Caluclate Functional Connectivity
% returns Functional Connectivity (FC) and p-values (P)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [FC, P] = calcFunctionalConnectivity(X, exSignal, nodeControl, exControl, isFullNode)
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
    [FC, P] = corr(X.', X.');

    if ~isempty(exSignal)
        FC = FC(1:nodeNum,:); P = P(1:nodeNum,:); 
    end
    if ~isempty(nodeControl)
        nodeControl(nodeControl==0) = nan;
        FC(:,1:nodeNum) = FC(:,1:nodeNum) .* nodeControl;
        P(:,1:nodeNum) = P(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl(exControl==0) = nan;
        FC(:,nodeNum+1:end) = FC(:,nodeNum+1:end) .* exControl;
        P(:,nodeNum+1:end) = P(:,nodeNum+1:end) .* exControl;
    end
end
