%%
% Plot Support Vector Normalized Partial Cross-Correlation (SvNPCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  kernel       kernel for SVM (default:'linear', 'gaussian', 'rbf')
%  kernelScale  kernelScale for SVM (default:'auto', 1)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [NPCC, lags] = plotSvPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, kernel, kernelScale, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, kernelScale = 'auto'; end
    if nargin < 6, kernel = 'linear'; end
    if nargin < 5, maxlag = 5; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    
    [NPCC, lags] = calcSvPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, kernel, kernelScale, isFullNode);

    % show coherence matrix
    fnum = size(NPCC,3);
    n = ceil(sqrt(fnum));
    for i=1:fnum
        subplot(n,n,i)
        clims = [-1,1];
        imagesc(NPCC(:,:,i),clims);
        daspect([1 1 1]);
        title(['SvPartial XCorr (' num2str(i-maxlag) ')']);
        xlabel('Source Nodes');
        ylabel('Target Nodes');
    end
    colorbar;
end
