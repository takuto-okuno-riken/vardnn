%%
% Plot DirectLiNGAM
% returns estimated causality (Aest) and so on.
% input:
%  X       multivariate time series matrix (node x time series)
%  rowcut  cut bottom rows of result gcI matris (default:0)

% Before using this function, download Dlingam-1.2 codes from
% https://sites.google.com/site/sshimizu06/Dlingamcode
% and add a path "Dlingam-1.2" and sub folders. And also download kernel-ICA 1.2 code from
% https://www.di.ens.fr/~fbach/kernel-ica/index.htm
% and add a path "kernel-ica1_2" and sub folders.

function [Aest, Best, stdeest, ciest, kest] = plotDirectLiNGAM(X, rowcut)
    if nargin < 2
        rowcut = 0;
    end
    % show Wavelet Coherence
    [Aest, Best, stdeest, ciest, kest] = calcDirectLiNGAM(X);
    if rowcut>0, mWCS(end-rowcut+1:end,:) = []; end
    clims = [-1,1];
    imagesc(Aest,clims);
    daspect([1 1 1]);
    title('DirectLiNGAM');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
