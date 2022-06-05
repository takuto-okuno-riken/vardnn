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
%  usegpu       use gpu calculation (default:false)

function [NCC, lags] = calcCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, isFullNode, usegpu)
    if nargin < 7, usegpu = false; end
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
    if usegpu
        Y = gpuArray(single(Y)); % 'half' does not support
    end

    % check all same value or not
    for i=1:nodeMax
        Ulen(i) = length(unique(single(Y(i,:)))); % 'half' does not support
    end

    NCC = nan(nodeNum,nodeMax,maxlag*2+1,class(X));
%    for i=1:nodeNum
    parfor i=1:nodeNum
        A = nan(nodeMax,maxlag*2+1,class(X));
        if usegpu
            A = gpuArray(single(A)); % 'half' does not support
        end
        x = Y(i,:).';
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            if Ulen(i)==1 && Ulen(j)==1
                A(j,:) = 0;
            else
                [A(j,:), ~] = xcov(x,Y(j,:).',maxlag,'normalized');
            end
        end
        if usegpu
            NCC(i,:,:) = gather(A);
        else
            NCC(i,:,:) = A;
        end
    end
    A = ones(nodeNum,'logical'); A = tril(A,-1);
    idx = find(A==1);
    for i=1:size(NCC,3)
        B = NCC(:,1:nodeNum,i); C = B';
        B(idx) = C(idx);
        NCC(:,1:nodeNum,i) = B;
    end
    lags = -maxlag:maxlag;

    % output control
    NCC = NCC(1:nodeNum,:,:);
    if isFullNode == 0
        NCC = NCC(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        NCC(:,1:nodeNum,:) = NCC(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        NCC(:,nodeNum+1:end,:) = NCC(:,nodeNum+1:end,:) .* exControl;
    end
end
