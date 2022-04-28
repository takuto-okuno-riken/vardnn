%%
% Plot Convergent Cross Mapping Pairwise Granger Causality
% returns CCM causality (CCM) and p-values (P)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  E            embedding dimension (default:3)
%  tau          time delay used in the phase-space reconstruction (default:1)
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [CCM, P] = plotConvCrossMapPGC(X, exSignal, nodeControl, exControl, E, tau, range, alpha, isFullNode)
    if nargin < 9, isFullNode = 0; end
    if nargin < 8, alpha = 0.05; end
    if nargin < 7, range = 10; end
    if nargin < 6, tau = 1; end
    if nargin < 5, E = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [CCM, P] = calcConvCrossMapPGC_(X, exSignal, nodeControl, exControl, E, tau, alpha, isFullNode);
    if range <= 0
        sigma = std(CCM(:),1,'omitnan');
        avg = mean(CCM(:),'omitnan');
        gcI2 = (CCM - avg) / sigma;
        range = 3;
    else
        gcI2 = CCM;
    end
    clims = [-range, range];
    imagesc(gcI2,clims);
    daspect([1 1 1]);
    title(['Convergent Cross Mapping pairwise GC Index (E=' num2str(E) ', tau=' num2str(tau) ')']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
