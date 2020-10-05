% Before using this function, download lingam-1.4.2 codes from
% http://www.cs.helsinki.fi/group/neuroinf/lingam/lingam.tar.gz
% and add a path "lingam-1.4.2" and sub folders, then remove a path "lingam-1.4.2/FastICA_21_octave" folders

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

[B,stde,ci,k] = lingam(X);

% plot matrix
figure;
clims = [-1 1];
imagesc(B,clims);
title('LiNGAM');
colorbar;

