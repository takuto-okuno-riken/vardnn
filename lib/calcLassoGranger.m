%%
% Calculate LassoGranger
% returns estimated causality (cause) and so on.
% input:
%  X            multivariate time series matrix (node x time series)
%  lags         number of lags for autoregression (default:1)
%  lambda       number of lags for autoregression (default:0.01)

% Before using this function, download Granger-causality codes from
% https://github.com/USC-Melady/Granger-causality
% and add a path "Granger-causality-master" and sub folders. And also download glmnet_matlab code from
% https://github.com/growlix/glmnet_matlab
% and add a path "glmnet_matlab-master" and sub folders. (Original glmnet was for windows7 and mex do not work)

function [cause] = calcLassoGranger(X, lags, lambda)
    if nargin < 3, lambda = 1e-2; end
    if nargin < 2, lags = 1; end
    N = size(X,1);
    cause = zeros(N, N);
    for in = 1:N
        index = [in, 1:(in-1), (in+1):N];
        [~, temp] = lassoGranger(X(index, :), lags, lambda, 'l');
        cause(in, :) = temp([2:in, 1, (in+1):N])';
    end
end
