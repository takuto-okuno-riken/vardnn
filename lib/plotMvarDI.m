%%
% Plotting mVAR (multivaliate Vector Auto-Regression) DI matrix
% returns mVAR (multivaliate Vector Auto-Regression) DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          mVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  range        plotting minimum and maximum range of DI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (optional)

function [DI, DIsub, coeff] = plotMvarDI(net, nodeControl, exControl, range, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, range = 10; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    [DI, DIsub, coeff] = calcMvarDI(net, nodeControl, exControl, isFullNode);
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
    title('mVAR Directional Influence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
