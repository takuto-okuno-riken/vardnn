%%
% Caluclate normalized cross-correlation (faster version)
% returns normalized cross-correlation (NCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  isFullNode   return both node & exogenous causality matrix (optional)
%  usegpu       use gpu calculation (default:false)

function [NCC, lags] = calcCrossCorrelation_(X, exSignal, nodeControl, exControl, maxlag, isFullNode, usegpu)
    if nargin < 7, usegpu = false; end
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, maxlag = 5; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    if usegpu
        Y = gpuArray(single(Y)); % 'half' does not support
    end

    % calc cross correlation
    NCC = nan(nodeMax,nodeMax,maxlag*2+1,class(X));
    if usegpu
        NCC = gpuArray(single(NCC)); % 'half' does not support
    end
    [C, lags] = xcov(Y',maxlag,'normalized');
    parfor i=1:size(C,1)
        NCC(:,:,i) = reshape(C(i,:),nodeMax,nodeMax)';
    end
    
    if usegpu
        A = gather(NCC);
        NCC = nan(nodeMax,nodeMax,maxlag*2+1,class(X));
        NCC(:,:,:) = A(:,:,:);
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
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        NCC(:,nodeNum+1:end,:) = NCC(:,nodeNum+1:end,:) .* exControl;
    end
end
