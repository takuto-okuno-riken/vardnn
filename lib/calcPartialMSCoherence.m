%%
% Caluclate Partial Magnitude-Squared Coherence (MSC)
% returns Partial Magnitude-Squared Coherence matrix (MSC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  nfft         sampling number of DFT (default:20)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [MSC, f] = calcPartialMSCoherence(X, exSignal, nodeControl, exControl, nfft, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, nfft = 20; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    
    if mod(nfft,2)==0, n = nfft/2 + 1; else n = (nfft+1)/2; end
    
    MSC = nan(nodeNum,nodeMax,n);
    fullIdx = 1:nodeMax;
    for i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);
        
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
            x = Y(i,:).';
            y = Y(j,:).';
            idx = setdiff(nodeIdx,j);
            z = [Y(idx,:).', ones(sigLen,1)];

            [b1,bint1,r1] = regress(x, z);
            [b2,bint2,r2] = regress(y, z);

            [MSC(i,j,:), f] = mscohere(r1,r2,[],[],nfft);
            MSC(j,i,:) = MSC(i,j,:);
        end
    end

    % output control
    MSC = MSC(1:nodeNum,:,:);
    if isFullNode == 0
        MSC = MSC(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        MSC(:,1:nodeNum,:) = MSC(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        MSC(:,nodeNum+1:end,:) = MSC(:,nodeNum+1:end,:) .* exControl;
    end
end
