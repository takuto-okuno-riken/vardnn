%%
% convert BOLD signal (gaussian distribution) to neural network signal (uniform distribution) [0, 1]
% return Y is time series matrix (node x time series)
% input:
%  X          multivariate time series matrix (node x time series)

function [Y, sig, m, maxsi, minsi] = gaussian2uniformSignal(X)
    maxsi = max(X(:));
    minsi = min(X(:));
    sig = sqrt(var(X(:),1));
    m = mean(X(:));
    Xn = (X-m)/ sig;
    Y = 0.5 * erfc(-Xn / sqrt(2));
end
