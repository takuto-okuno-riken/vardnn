%%
% Caluclate Wavelet Coherence
% returns mean Wavelet Cross Spectol matrix (mWCS), full set of Wavelet
% Coherence (WCOH) and Cross Spectols (WCS)
% input:
%  X     multivariate time series matrix (node x time series)

function [mWCS, WCOH, WCS] = calcWaveletCoherence(X)
    n = size(X,1);
    mWCS = ones(n,n);
    WCOH = cell(n,n);
    WCS = cell(n,n);
    for i=1:n
        for j=1:n
            if i==j, continue; end
            [wcoh, wcs] = wcoherence(X(i,:),X(j,:));
            WCOH{i,j} = wcoh;
            WCS{i,j} = wcs;
            mWCS(i,j) = nanmean(real(wcs),'all');
        end
    end
end
