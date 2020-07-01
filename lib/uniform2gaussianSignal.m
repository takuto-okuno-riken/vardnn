%%
% convert neural network signal (uniform distribution) [0, 1] to BOLD signal (gaussian distribution)
% return X is time series matrix (node x time series)
% input:
%  Y          multivariate time series matrix (node x time series)
%  sig        sigma value of original X matrix
%  m          mean value of original X matrix
%  maxsi      maximum value of original X matrix
%  minsi      minimum value of original X matrix
function X = uniform2gaussianSignal(Y, sig, m, maxsi, minsi)
    y2 = -erfcinv(2 * Y) * sqrt(2);
    X = sig * y2 + m;
    % replace inf to original max or min values
    idx = find(isinf(X)&X>0);
    X(idx) = maxsi;
    idx = find(isinf(X)&X<0);
    X(idx) = minsi;
end
