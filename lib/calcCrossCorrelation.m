%%
% Caluclate normalized cross-correlation
% returns normalized cross-correlation (NCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [NCC, lags] = calcCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, isFullNode)
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

    NCC = nan(nodeNum,nodeMax,maxlag*2+1);
    for i=1:nodeNum
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            if Ulen(i)==1 && Ulen(j)==1
                NCC(i,j,:) = 0;
                lags = -maxlag:maxlag;
            else
                x = Y(i,:).';
                y = Y(j,:).';
                [NCC(i,j,:), lags] = xcov(x,y,maxlag,'normalized');
            end
            NCC(j,i,:) = NCC(i,j,:);
        end
    end

    % output control
    NCC = NCC(1:nodeNum,:,:);
    if isFullNode == 0
        NCC = NCC(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        NCC(:,1:nodeNum,:) = NCC(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        NCC(:,nodeNum+1:end,:) = NCC(:,nodeNum+1:end,:) .* exControl;
    end
end
