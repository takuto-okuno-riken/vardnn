%%
% Calculate cross-power spectal density
% returns cross-power spectal density matrix (CPM)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  nfft         sampling number of DFT (default:20)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [CPM] = calcCPSD(X, exSignal, nodeControl, exControl, nfft, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, nfft = 20; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    
    if mod(nfft,2)==0, n = nfft/2 + 1; else n = (nfft+1)/2; end
    
    CPM = nan(nodeNum,nodeMax,n);
    for i=1:nodeNum
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            x = Y(i,:).';
            y = Y(j,:).';
            CPM(i,j,:) = cpsd(x,y,[],[],nfft);
            CPM(j,i,:) = CPM(i,j,:);
        end
    end

    % output control
    CPM = CPM(1:nodeNum,:,:);
    if isFullNode == 0
        CPM = CPM(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        CPM(:,1:nodeNum,:) = CPM(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        CPM(:,nodeNum+1:end,:) = CPM(:,nodeNum+1:end,:) .* exControl;
    end
end
