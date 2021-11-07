%%
% Plot Random Forest Partial Correlation
% returns Random Forest Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  numTrees     number of trees for random forest (default:10)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = plotRfPartialCorrelation(X, exSignal, nodeControl, exControl, numTrees, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, numTrees = 10; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    PC = calcRfPartialCorrelation(X, exSignal, nodeControl, exControl, numTrees, isFullNode);

    % show partial correlation
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title(['Random Forest Partial Correlation (numTrees=' num2str(numTrees) ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
