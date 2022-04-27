%%
% Plot Convergent Cross Mapping - FC (subtract FC)
% returns CCM causality (CCM), FC p-values (Pfc) and CCM p-values (Pccm).
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  E            embedding dimension (default:3)
%  tau          time delay used in the phase-space reconstruction (default:1)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [CCM, Pfc, Pccm] = plotConvCrossMapSubFC(X, exSignal, nodeControl, exControl, E, tau, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, tau = 1; end
    if nargin < 5, E = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [CCM, Pfc, Pccm] = calcConvCrossMapSubFC(X, exSignal, nodeControl, exControl, E, tau, isFullNode);
    clims = [-1,1];
    imagesc(CCM,clims);
    daspect([1 1 1]);
    title(['Convergent Cross Mapping - FC (E=' num2str(E) ', tau=' num2str(tau) ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
