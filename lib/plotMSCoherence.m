%%
% plot Magnitude-Squared Coherence (MSC) matrix
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  nfft         sampling number of DFT (default:20)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [MSC, f] = plotMSCoherence(X, exSignal, nodeControl, exControl, nfft, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, nfft = 20; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    
    [MSC, f] = calcMSCoherence(X, exSignal, nodeControl, exControl, nfft, isFullNode);

    % show coherence matrix
    fnum = size(MSC,3);
    n = ceil(sqrt(fnum));
    for i=1:fnum
        subplot(n,n,i)
        clims = [0,1];
        imagesc(MSC(:,:,i),clims);
        daspect([1 1 1]);
        title(['MS Coherence (' num2str(i) ')']);
        xlabel('Source Nodes');
        ylabel('Target Nodes');
    end
    colorbar;
end
