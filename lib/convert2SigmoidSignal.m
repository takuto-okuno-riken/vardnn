%%
% convert signal. normalized by sigma, then sigmoid [0, 1]
% return Y is time series matrix (node x time series)
% input:
%  X          multivariate time series matrix (node x time series)
%  a          sigmoidal 'a' coefficient (default:1)

function [Y, sig, m, maxsi, minsi] = convert2SigmoidSignal(X, a)
    if nargin < 2
        a = 1;
    end
    maxsi = max(X(:));
    minsi = min(X(:));
    sig = sqrt(var(X(:)));
    m = mean(X(:));
    Xn = (X-m)/ sig;
    Y = sigmf(Xn,[a 0]);
end
