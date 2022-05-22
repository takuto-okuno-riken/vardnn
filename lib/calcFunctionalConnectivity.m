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
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % set node input
    X = [X; exSignal];

    % check all same value or not
    for i=1:nodeNum+exNum
        A = unique(X(i,:));
        if length(A)==1
            X(i,mod(i,sigLen)+1) = A + 1.0e-8;
        end
    end
    
    [FC, P] = corr(X.', X.');

    % output control
    FC = FC(1:nodeNum,:); P = P(1:nodeNum,:); 
    if isFullNode == 0
        FC = FC(:,1:nodeNum); P = P(:,1:nodeNum); 
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        FC(:,1:nodeNum) = FC(:,1:nodeNum) .* nodeControl;
        P(:,1:nodeNum) = P(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        FC(:,nodeNum+1:end) = FC(:,nodeNum+1:end) .* exControl;
        P(:,nodeNum+1:end) = P(:,nodeNum+1:end) .* exControl;
    end
end
