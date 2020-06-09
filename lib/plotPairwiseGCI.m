%%
% Plotting pairwise Granger causality Index matrix
% returns Granger Causality Index (gcI)
% input:
%  X       multivariate time series matrix (node x time series)
%  lags    number of lags for autoregression (default:3)
%  range   plotting minimum and maximum range of GCI (default:10)
%  rowcut  cut bottom rows of result gcI matris (default:0)

function [gcI] = plotPairwiseGCI(X, lag, range, rowcut)
    if nargin < 4
        rowcut = 0;
    end
    if nargin < 3
        range = 10;
    end
    if nargin < 2
        lag = 3;
    end
    gcI = calcPairwiseGCI(X, lag);
    if rowcut>0, gcI(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(gcI,clims);
    daspect([1 1 1]);
    title('Granger Causality Index (GC-EC)');
    colorbar;
end
