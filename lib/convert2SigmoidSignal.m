%%
% convert signal. normalized by sigma, then sigmoid [0, 1]
% return Y is time series matrix (node x time series)
% input:
%  X          multivariate time series matrix (node x time series)
%  centroid   signal centroid value for normalization (i.e. BOLD: set 0, other:auto)(option)
%  a          sigmoidal 'a' coefficient (default:1)

function [Y, sig, c, maxsi, minsi] = convert2SigmoidSignal(X, centroid, a)
    if nargin < 3
        a = 1;
    end
    if nargin < 2
        centroid = NaN;
    end
    maxsi = max(X(:));
    minsi = min(X(:));
    sig = sqrt(var(X(:)));
    if isnan(centroid)
        c = mean(X(:));
    else
        c = centroid;
    end
    Xn = (X-c)/ sig;
    Y = sigmf(Xn,[a 0]);
end
