%%
% Plotting Pairwised VAR DNN-DI matrix
% returns Pairwised VAR DNN DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained Pairwised VAR DNN network structure
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  range        plotting minimum and maximum range of DI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (optional)

function [DI, DIsub] = plotPvarDnnDI(net, nodeControl, exControl, range, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, range = 10; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    [DI, DIsub] = calcPvarDnnDI(net, nodeControl, exControl, isFullNode);
    if range <= 0
        sigma = std(DI(:),1,'omitnan');
        avg = mean(DI(:),'omitnan');
        DI2 = (DI - avg) / sigma;
        range = 3;
    else
        DI2 = DI;
    end
    clims = [-range, range];
    imagesc(DI2,clims);
    daspect([1 1 1]);
    title('pairwise VAR DNN Directional Influence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
