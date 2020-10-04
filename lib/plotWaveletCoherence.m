%%
% Plot Wavelet Coherence
% returns mean Wavelet Cross Spectol matrix (mWCS), full set of Wavelet
% Coherence (WCOH) and Cross Spectols (WCS)
% input:
%  X       multivariate time series matrix (node x time series)
%  rowcut  cut bottom rows of result gcI matris (default:0)

function [mWCS, WCOH, WCS] = plotWaveletCoherence(X, rowcut)
    if nargin < 2
        rowcut = 0;
    end
    % show Wavelet Coherence
    [mWCS, WCOH, WCS] = calcWaveletCoherence(X);
    if rowcut>0, mWCS(end-rowcut+1:end,:) = []; end
    clims = [-1,1];
    imagesc(mWCS,clims);
    daspect([1 1 1]);
    title('Wavelet Coherence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
