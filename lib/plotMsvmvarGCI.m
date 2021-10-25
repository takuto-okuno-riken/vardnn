%%
% Plotting mSvmVAR (multivaliate SVM Vector Auto-Regression) Granger Causality
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mSvmVAR network
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotMsvmvarGCI(X, exSignal, nodeControl, exControl, net, range, alpha, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, alpha = 0.05; end
    if nargin < 6, range = 10; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMsvmvarGCI(X, exSignal, nodeControl, exControl, net, alpha, isFullNode);
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
    title('mSvmVAR Granger Causality Index');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
