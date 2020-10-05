%%
% convert signal. normalized by sigma, then sigmoid [0, 1]
% return Y is time series matrix (node x time series)
% input:
%  X          multivariate time series matrix (node x time series)
%  signalm    signal mean value for normalization (i.e. BOLD: set 0, other:auto)(option)
%  a          sigmoidal 'a' coefficient (default:1)

function [Y, sig, m, maxsi, minsi] = convert2SigmoidSignal(X, signalm, a)
    if nargin < 3
        a = 1;
    end
    if nargin < 2
        signalm = NaN;
    end
    maxsi = max(X(:));
    minsi = min(X(:));
    sig = sqrt(var(X(:)));
    if isnan(signalm)
        m = mean(X(:));
    else
        m = signalm;
    end
    Xn = (X-m)/ sig;
    Y = sigmf(Xn,[a 0]);
end
