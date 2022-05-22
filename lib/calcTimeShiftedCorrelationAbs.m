%%
% Caluclate time shifted Correlation (Abs) (Functional Connectivity)
% cross-correlation type signal processing
% returns Functional Connectivity (FC) and p-values (P)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for time shift (default:3)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [FC, P] = calcTimeShiftedCorrelationAbs(X, exSignal, nodeControl, exControl, lags, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    m = size(X,2);

    X = [X; exSignal];
    
    if lags > 0
        X1 = repmat(X(:,1:m-lags), [1 lags]);
        X2 = [];
        for i=1:lags
            X2 = [X2, X(:,1+i:m-(lags-i))];
        end
        [FC, P] = corr(X2.', X1.');
    else
        [FC, P] = corr(X.', X.');
    end
    FC = abs(FC);

    % output control
    if ~isempty(exSignal)
        FC = FC(1:nodeNum,:); P = P(1:nodeNum,:); 
    end
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
