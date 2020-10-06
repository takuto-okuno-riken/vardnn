%%
% convert sigmoidal signal [0, 1] and denormalized to original signal 
% return X is time series matrix (node x time series)
% input:
%  Y          multivariate time series matrix (node x time series)
%  sig        sigma value of original X matrix
%  c          centroid value of original X matrix
%  maxsi      maximum value of original X matrix
%  minsi      minimum value of original X matrix
%  a          sigmoidal 'a' coefficient (default:1)

function X = convert2InvSigmoidSignal(Y, sig, c, maxsi, minsi, a)
    if nargin < 6
        a = 1;
    end
    y2 = real(log(Y ./ (Y-1)) / a);
    X = sig * y2 + c;
    % replace inf to original max or min values
    idx = find(isinf(X)&X>0);
    X(idx) = maxsi;
    idx = find(isinf(X)&X<0);
    X(idx) = minsi;
end
