%%
% Plotting DLCM Granger causality Index matrix
% returns DLCM Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotDlcmGCI(X, exSignal, nodeControl, exControl, netDLCM, range, alpha, isFullNode)
    if nargin < 8
        isFullNode = 0;
    end
    if nargin < 7
        alpha = 0.05;
    end
    if nargin < 6
        range = 10;
    end
    [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI(X, exSignal, nodeControl, exControl, netDLCM, alpha, isFullNode);
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
    title('DLCM Granger Causality');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
