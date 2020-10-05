% Before using this function, download Dlingam-1.2 codes from
% https://sites.google.com/site/sshimizu06/Dlingamcode
% and add a path "Dlingam-1.2" and sub folders. And also download kernel-ICA 1.2 code from
% https://www.di.ens.fr/~fbach/kernel-ica/index.htm
% and add a path "kernel-ica1_2" and sub folders.

nodeNum = 8;  % node number
sigLen = 100; % signal length

% generate random signals
E = randn(nodeNum, sigLen); 

% synchronize signal 6 -> 2, 6 invert 3, 7 invert 4
B = zeros(nodeNum, nodeNum);
B(2,6) = 1;
B(3,6) = -1; B(4,7) = -1;
A = inv( eye(nodeNum) - B );
X = A * E;

[Best, stdeest, ciest, kest] = Dlingam(X);
Aest = estA(X,kest); % or inv(eye(p)-Best);

% plot matrix
figure;
clims = [-1 1];
imagesc(Aest,clims);
title('DirectLiNGAM');
colorbar;

