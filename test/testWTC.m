nodeNum = 8;  % node number
sigLen = 100; % signal length

% generate random signals
X = rand(nodeNum, sigLen); 

% synchronize signal 6 == 2, 6 invert 3, 7 invert 4
X(2,:) = X(6,:);
X(3,:) = 1 - X(6,:);
X(4,:) = 1 - X(7,:);

[wcoh1,wcs1] = wcoherence(X(2,:),X(6,:));

[wcoh2,wcs2] = wcoherence(X(2,:),X(3,:));

[mWCS, WCOH, WCS] = calcWaveletCoherence(X);

% plot matrix
figure;
clims = [-1 1];
imagesc(mWCS,clims);
title('mean Wavelet Coherence');
colorbar;

