%%
% Calculate Discrete Fourier transform
% returns One-sided amplitude spectrum D (node x sampling spectrum)
% input:
%  X         multivariate time series matrix (node x time series)
%  n         DFT sampling number (even number) (default: 100)
%  Fs        sampling frequency of time seriese (default: 0.5)
function [D, P1] = calcDft(X, n, Fs)
    if nargin < 3, Fs = 0.5; end
    if nargin < 2, n = 100; end

    Y = fft(X,n,2);
    P2 = abs(Y/n);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    D = P1(:,2:end-1);
end
