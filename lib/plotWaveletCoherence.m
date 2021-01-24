%%
% Plot Wavelet Coherence
% returns mean Wavelet Cross Spectol matrix (mWCS), full set of Wavelet
% Coherence (WCOH) and Cross Spectols (WCS)
% input:
%  X       multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [mWCS, WCOH, WCS] = plotWaveletCoherence(X, exSignal, nodeControl, exControl, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [mWCS, WCOH, WCS] = calcWaveletCoherence(X, exSignal, nodeControl, exControl, isFullNode);

    % show Wavelet Coherence
    clims = [-1,1];
    imagesc(mWCS,clims);
    daspect([1 1 1]);
    title('Wavelet Coherence');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
