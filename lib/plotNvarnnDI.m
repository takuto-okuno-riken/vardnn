%%
% plot nVARNN Directional Influence matrix
% input:
%  net          trained nVARNN network
%  nodeControl  node control matrix (1 x node) (optional)
%  exControl    exogenous input control matrix for each node (1 x exogenous input) (optional)
%  range        plotting minimum and maximum range of DI (default:0.5)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = plotNvarnnDI(net, nodeControl, exControl, range, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, range = 0.5; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    [DI, DIsub] = calcNvarnnDI(net, nodeControl, exControl, isFullNode);
    % show nVARNN weight causality of predicted node signals
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
    title('nVARNN Directional Influence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
