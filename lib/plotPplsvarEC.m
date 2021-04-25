%%
% Plotting pPLSVAR (pairwise PLS Vector Auto-Regression) EC matrix
% returns pPLSVAR (pairwise PLS Vector Auto-Regression) EC matrix (EC) and impaired node signals (ECsub)
% input:
%  net          pPLSVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (optional)

function [EC, ECsub] = plotPplsvarEC(net, nodeControl, exControl, range, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, range = 10; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    [EC, ECsub] = calcPplsvarEC(net, nodeControl, exControl, isFullNode);
    if range <= 0
        sigma = std(EC(:),1,'omitnan');
        avg = mean(EC(:),'omitnan');
        EC2 = (EC - avg) / sigma;
        range = 3;
    else
        EC2 = EC;
    end
    clims = [-range, range];
    imagesc(EC2,clims);
    daspect([1 1 1]);
    title('pPLSVAR Effective Connectivity');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
