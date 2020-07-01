%%
% convert BOLD signal to neural network input signal [0, 1]
% return Y is time series matrix (node x time series)
% input:
%  X          multivariate time series matrix (node x time series)
%  a          sigmoidal 'a' coefficient (default:4)

function Y = bold2dnnSignal(X, a)
    if nargin < 2
        a = 4;
    end
    Y = sigmf(X,[a 0]);
end
