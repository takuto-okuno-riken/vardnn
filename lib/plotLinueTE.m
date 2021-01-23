%%
% Plotting LINER-Uniform Embedding Transfer Entropy (LINUE-TE) matrix
% returns Transfer Entropy matrix (TE), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
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

function [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotLinueTE(X, exSignal, nodeControl, exControl, lags, range, alpha, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, alpha = 0.05; end
    if nargin < 6, range = 10; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcLinueTE(X, exSignal, nodeControl, exControl, lags, alpha, isFullNode);
    clims = [0, range];
    if range <= 0
        sigma = std(TE(:),1,'omitnan');
        avg = mean(TE(:),'omitnan');
        TE2 = (TE - avg) / sigma;
        clims = [-3, 3];
    else
        TE2 = TE;
    end
    imagesc(TE2,clims);
    daspect([1 1 1]);
    title('Transfer Entropy (LINUE)');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
