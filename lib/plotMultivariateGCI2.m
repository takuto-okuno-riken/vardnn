%%
% Plotting multivariate Granger causality Index matrix
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X       multivariate time series matrix (node x time series)
%  lags    number of lags for autoregression (default:3)
%  range   plotting minimum and maximum range of GCI (default:10)
%          if range==0, range shows standard deviation [-5 sigma, 5 sigma]
%  rowcut  cut bottom rows of result gcI matris (default:0)
%  alpha  the significance level of F-statistic (default:0.05)

function [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotMultivariateGCI2(X, lag, range, rowcut, alpha)
    if nargin < 5
        alpha = 0.05;
    end
    if nargin < 4
        rowcut = 0;
    end
    if nargin < 3
        range = 10;
    end
    if nargin < 2
        lag = 3;
    end
   [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMultivariateGCI2(X, lag, alpha);
    if range <= 0
        sigma = std(gcI(:),'omitnan');
        avg = mean(gcI(:),'omitnan');
        gcI = (gcI - avg) / sigma;
        range = 5;
    end
    if rowcut>0, gcI(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(gcI,clims);
    daspect([1 1 1]);
    title('multivariate Granger Causality Index');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
