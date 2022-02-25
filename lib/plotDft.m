%%
% Plot Discrete Fourier transform
% returns One-sided amplitude spectrum D (node x sampling spectrum)
% input:
%  X         multivariate time series matrix (node x time series)
%  n         DFT sampling number (even number) (default: 100)
%  Fs        sampling frequency of time seriese (default: 0.5)

function [D, P1] = plotDft(X, n, Fs)
    if nargin < 3, Fs = 0.5; end
    if nargin < 2, n = 100; end

    [D, P1] = calcDft(X, n, Fs);

    % plot Amplitude Spectrum
    f = Fs*(0:(n/2))/n;
    plot(f, P1.');
    title('Single-Sided Amplitude Spectrum of S(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
end
