%%
% Plot Principal Component Partial Correlation
% returns Principal Component Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  explainedTh  explained threshold for PCA components (default:0.99)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = plotPcPartialCorrelation(X, exSignal, nodeControl, exControl, explainedTh, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, explainedTh = 0.99; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    PC = calcPcPartialCorrelation(X, exSignal, nodeControl, exControl, explainedTh, isFullNode);

    % show partial correlation
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title(['PC Partial Correlation (th=' num2str(explainedTh) ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
