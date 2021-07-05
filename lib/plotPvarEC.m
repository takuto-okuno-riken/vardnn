%%
% Plotting PVAR (pairwised vector auto-regression) EC matrix
% returns PVAR (pairwised vector auto-regression) EC matrix (EC) and impaired node signals (ECsub)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [EC, ECsub, coeff] = plotPvarEC(X, exSignal, nodeControl, exControl, lags, range, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, range = 10; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [EC, ECsub, coeff] = calcPvarEC(X, exSignal, nodeControl, exControl, lags, isFullNode);
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
    title('PVAR Directional Influence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
