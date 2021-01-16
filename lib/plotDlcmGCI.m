%%
% Plotting DLCM Granger causality Index matrix
% returns DLCM Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  inSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network
%  range        plotting minimum and maximum range of GCI (default:10)
%  rowcut       cut bottom rows of result gcI matris (default:0)
%  alpha        the significance level of F-statistic (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM, range, rowcut, alpha, isFullNode)
    if nargin < 9
        isFullNode = 0;
    end
    if nargin < 8
        alpha = 0.05;
    end
    if nargin < 7
        rowcut = 0;
    end
    if nargin < 6
        range = 10;
    end
    [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM, alpha, isFullNode);
    if range <= 0
        sigma = std(gcI(:),1,'omitnan');
        avg = mean(gcI(:),'omitnan');
        gcI = (gcI - avg) / sigma;
        range = 3;
    end
    if rowcut>0, gcI(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(gcI,clims);
    daspect([1 1 1]);
    title('DLCM Granger Causality');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
