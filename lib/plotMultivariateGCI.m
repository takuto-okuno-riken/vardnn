%%
% Plotting multivariate Granger causality Index matrix
% returns Granger Causality Index (gcI)
% input:
%  X       multivariate time series matrix (node x time series)
%  lags    number of lags for autoregression (default:3)
%  range   plotting minimum and maximum range of GCI (default:10)
%          if range==0, range shows standard deviation [-5 sigma, 5 sigma]
%  rowcut  cut bottom rows of result gcI matris (default:0)

function [gcI] = plotMultivariateGCI(X, lag, range, rowcut)
    if nargin < 4
        rowcut = 0;
    end
    if nargin < 3
        range = 10;
    end
    if nargin < 2
        lag = 3;
    end
    gcI = calcMultivariateGCI(X, lag);
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
    title('multivariate Granger Causality Index');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
