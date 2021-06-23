%%
% Plotting PLassoVAR (pairwised Lasso vector auto-regression) EC matrix
% returns PLassoVAR (pairwised Lasso vector auto-regression) EC matrix (EC) and impaired node signals (ECsub)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  lambda       lambda for the Lasso (default:0.01)
%  elaAlpha     Elastic Net Alpha for the Lasso (default:1)
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [EC, ECsub, coeff] = plotPlassovarEC(X, exSignal, nodeControl, exControl, lags, lambda, elaAlpha, range, isFullNode)
    if nargin < 9, isFullNode = 0; end
    if nargin < 8, range = 10; end
    if nargin < 7, elaAlpha = 1; end
    if nargin < 6, lambda = 0.01; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [EC, ECsub, coeff] = calcPlassovarEC(X, exSignal, nodeControl, exControl, lags, lambda, elaAlpha, isFullNode);
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
    title('PLassoVAR Effective Connectivity');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
