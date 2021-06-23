%%
% Plot LassoGranger
% returns estimated causality (cause) and so on.
% input:
%  X            multivariate time series matrix (node x time series)
%  lags         number of lags for autoregression (default:1)
%  lambda       number of lags for autoregression (default:0.01)
%  range        plotting minimum and maximum range of GCI (default:1)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]

% Before using this function, download Granger-causality codes from
% https://github.com/USC-Melady/Granger-causality
% and add a path "Granger-causality-master" and sub folders. And also download glmnet_matlab code from
% https://github.com/growlix/glmnet_matlab
% and add a path "glmnet_matlab-master" and sub folders. (Original glmnet was for windows7 and mex do not work)

function [EC] = plotLassoGranger(X, lags, lambda, range)
    if nargin < 4, range = 1; end
    if nargin < 3, lambda = 1e-2; end
    if nargin < 2, lags = 1; end
    % show LassoGranger
    EC = calcLassoGranger(X, lags, lambda);
    if range <= 0
        sigma = std(EC(:),1,'omitnan');
        avg = mean(EC(:),'omitnan');
        EC2 = (EC - avg) / sigma;
        range = 3;
    else
        EC2 = EC;
    end
    clims = [-range,range];
    imagesc(EC2,clims);
    daspect([1 1 1]);
    title('LassoGranger');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
