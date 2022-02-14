%%
% Plot Normalized Partial Cross-Correlation (NPCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [NPCC, lags] = plotPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, maxlag = 5; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    
    [NPCC, lags] = calcPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, isFullNode);

    % show coherence matrix
    fnum = size(NPCC,3);
    n = ceil(sqrt(fnum));
    for i=1:fnum
        subplot(n,n,i)
        clims = [-1,1];
        imagesc(NPCC(:,:,i),clims);
        daspect([1 1 1]);
        title(['Partial XCorr (' num2str(i-maxlag) ')']);
        xlabel('Source Nodes');
        ylabel('Target Nodes');
    end
    colorbar;
end
