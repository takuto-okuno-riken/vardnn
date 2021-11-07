%%
% Plot DNN Partial Correlation
% returns DNN Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  options         training options
%  activateFunc    activation function for each layer (default:@reluLayer)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = plotDnnPartialCorrelation(X, exSignal, nodeControl, exControl, options, activateFunc, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, activateFunc = @reluLayer; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    PC = calcDnnPartialCorrelation(X, exSignal, nodeControl, exControl, options, activateFunc, isFullNode);

    % show partial correlation
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title(['DNN Partial Correlation']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
