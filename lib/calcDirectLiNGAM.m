%%
% Calculate DirectLiNGAM
% returns estimated causality (Aest) and so on.
% input:
%  X     multivariate time series matrix (node x time series)

% Before using this function, download Dlingam-1.2 codes from
% https://sites.google.com/site/sshimizu06/Dlingamcode
% and add a path "Dlingam-1.2" and sub folders. And also download kernel-ICA 1.2 code from
% https://www.di.ens.fr/~fbach/kernel-ica/index.htm
% and add a path "kernel-ica1_2" and sub folders.

function [Aest, Best, stdeest, ciest, kest] = calcDirectLiNGAM(X)
    [Best, stdeest, ciest, kest] = Dlingam(X);
    Aest = estA(X,kest); % or inv(eye(p)-Best);
end
