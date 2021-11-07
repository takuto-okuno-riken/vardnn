%%
% Plot Gaussian Processes Partial Correlation
% returns Gaussian Processes Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  kernel          kernel for GP (default:'squaredexponential', 'ardsquaredexponential', ... please see https://jp.mathworks.com/help/stats/fitrgp.html )
%  basis           basis for GP (default:'linear', 'none', 'constant')
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = plotGpPartialCorrelation(X, exSignal, nodeControl, exControl, kernel, basis, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, basis = 'linear'; end
    if nargin < 5, kernel = 'squaredexponential'; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    PC = calcGpPartialCorrelation(X, exSignal, nodeControl, exControl, kernel, basis, isFullNode);

    % show partial correlation
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title(['GP Partial Correlation (kernel=' kernel ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
