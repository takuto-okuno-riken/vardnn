%%
% convert DCM y (BOLD) signal to neural network input signal [0, 1]
% return Y is time series matrix (node x time series)
% input:
%  X          multivariate time series matrix (node x time series)

function Y = dcmY2DlcmSignal(X)
    minSi = min(min(X));
    maxSi = max(max(X));
    maxSi = floor((maxSi-minSi)*10) / 10;
    minSi = ceil(minSi*10) / 10; % should be minus
    Y = (X - minSi) / maxSi;
end
