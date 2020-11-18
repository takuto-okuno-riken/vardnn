%%
% Plotting LINER-Uniform Embedding Transfer Entropy (LINUE-TE) matrix
% returns Transfer Entropy matrix (TE), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X       multivariate time series matrix (node x time series)
%  lags    number of lags for autoregression (default:3)
%  range   plotting minimum and maximum range of TE (default:0.1)
%          if range==0, range shows standard deviation [-5 sigma, 5 sigma]
%  rowcut  cut bottom rows of result gcI matris (default:0)
%  alpha   the significance level of F-statistic (default:0.05)

function [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotLinueTE(X, lag, range, rowcut, alpha)
    if nargin < 5
        alpha = 0.05;
    end
    if nargin < 4
        rowcut = 0;
    end
    if nargin < 3
        range = 0.1;
    end
    if nargin < 2
        lag = 3;
    end
    clims = [0, range];
    [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcLinueTE(X, lag, alpha);
    if range <= 0
        sigma = std(TE(:),1,'omitnan');
        avg = mean(TE(:),'omitnan');
        TE = (TE - avg) / sigma;
        clims = [-3, 3];
    end
    if rowcut>0, TE(end-rowcut+1:end,:) = []; end
    imagesc(TE,clims);
    daspect([1 1 1]);
    title('Transfer Entropy (LINUE)');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
