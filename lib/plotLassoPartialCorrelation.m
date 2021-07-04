%%
% Plot Lasso Partial Correlation
% returns Lasso Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lambda       lambda for the Lasso (default:0.01)
%  elaAlpha     Elastic Net Alpha for the Lasso (default:1)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = plotLassoPartialCorrelation(X, exSignal, nodeControl, exControl, lambda, elaAlpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, elaAlpha = 1; end
    if nargin < 5, lambda = 0.01; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    PC = calcLassoPartialCorrelation(X, exSignal, nodeControl, exControl, lambda, elaAlpha, isFullNode);

    % show partial correlation
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title(['Lasso Partial Correlation (r=' num2str(lambda) ', a=' num2str(elaAlpha) ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
