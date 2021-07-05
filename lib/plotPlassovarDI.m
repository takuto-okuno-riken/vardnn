%%
% Plotting PLassoVAR (pairwised Lasso vector auto-regression) DI matrix
% returns PLassoVAR (pairwised Lasso vector auto-regression) DI matrix (DI) and impaired node signals (DIsub)
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

function [DI, DIsub, coeff] = plotPlassovarDI(X, exSignal, nodeControl, exControl, lags, lambda, elaAlpha, range, isFullNode)
    if nargin < 9, isFullNode = 0; end
    if nargin < 8, range = 10; end
    if nargin < 7, elaAlpha = 1; end
    if nargin < 6, lambda = 0.01; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [DI, DIsub, coeff] = calcPlassovarDI(X, exSignal, nodeControl, exControl, lags, lambda, elaAlpha, isFullNode);
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
    title('PLassoVAR Directional Influence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
