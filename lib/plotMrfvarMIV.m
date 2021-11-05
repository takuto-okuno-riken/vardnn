%%
% Plotting mRFVAR Mean Impact Value (MIV) matrix
% returns mRFVAR MIV matrix (MIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (1 x node) (optional)
%  exControl    exogenous input control matrix for each node (1 x exogenous input) (optional)
%  net          trained mRFVAR network
%  range        plotting minimum and maximum range of DI (default:1)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function MIV = plotMrfvarMIV(X, exSignal, nodeControl, exControl, net, range, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, range = 1; end

    MIV = calcMrfvarMIV(X, exSignal, nodeControl, exControl, net, isFullNode);
    if range <= 0
        sigma = std(MIV(:),1,'omitnan');
        avg = mean(MIV(:),'omitnan');
        MIV2 = (MIV - avg) / sigma;
        range = 3;
    else
        MIV2 = MIV;
    end
    clims = [-range, range];
    imagesc(MIV2,clims);
    daspect([1 1 1]);
    title('mRFVAR Mean Impact Value');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
