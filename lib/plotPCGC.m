%%
% Plotting mPCVAR (multivaliate Principal Component Vector Auto-Regression) Granger Causality
% returns Granger causality index matrix (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F), the critical value from the F-distribution (cvFd)
% and AIC, BIC (of node vector)
% input:
%  X            multivariate time series matrix (node x time series)
%  lags         number of lags for autoregression (default:1)
%  ndRate       (of ndmax) as the value selected as the knee of the curve (default:0.8)
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]

% Before using this function, download PartiallyConditionedGrangerCausality codes from
% https://github.com/danielemarinazzo/PartiallyConditionedGrangerCausality
% and add a path "PartiallyConditionedGrangerCausality-master" and sub folders. 

function gcI2 = plotPCGC(X, lags, ndRate, range)
    if nargin < 4, range = 10; end
    if nargin < 3, ndRate = 0.8; end
    if nargin < 2, lags = 1; end

    gcI = calcPCGC(X, lags, ndRate);
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
    title('PCGC Index');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
