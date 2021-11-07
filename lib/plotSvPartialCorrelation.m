%%
% Plot Support Vector Partial Correlation
% returns Support Vector Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  kernel          kernel for SVM (default:'linear', 'gaussian', 'rbf')
%  kernelScale     kernelScale for SVM (default:'auto', 1)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = plotSvPartialCorrelation(X, exSignal, nodeControl, exControl, kernel, kernelScale, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, kernelScale = 'auto'; end
    if nargin < 5, kernel = 'linear'; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    PC = calcSvPartialCorrelation(X, exSignal, nodeControl, exControl, kernel, kernelScale, isFullNode);

    % show partial correlation
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title(['SV Partial Correlation (kernel=' kernel ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
