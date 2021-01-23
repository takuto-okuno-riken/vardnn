%%
% Plotting pairwise Granger causality Index matrix
% returns Granger causality index (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F) and the critical value from the F-distribution (cvFd)
% and AIC, BIC
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC] = plotPairwiseGCI(X, exSignal, nodeControl, exControl, lag, range, alpha, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, alpha = 0.05; end
    if nargin < 6, range = 10; end
    if nargin < 5, lag = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [gcI, h, P, F, cvFd, AIC, BIC] = calcPairwiseGCI(X, exSignal, nodeControl, exControl, lag, alpha, isFullNode);
    if range <= 0
        sigma = std(gcI(:),1,'omitnan');
        avg = mean(gcI(:),'omitnan');
        gcI2 = (gcI - avg) / sigma;
        range = 3;
    else
        gcI2 = gcI;
    end
    clims = [-range, range];
    imagesc(gcI2,clims);
    daspect([1 1 1]);
    title('pairwised Granger Causality Index');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
