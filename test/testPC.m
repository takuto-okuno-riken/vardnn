nodeNum = 8;  % node number
sigLen = 100; % signal length

% generate random signals
%X = rand(nodeNum, sigLen); 
load('test/testTrain-rand500-uniform.mat');
X = si(1:8, 1:sigLen);
x = X(1,:).';
y = X(2,:).';
z = X(3:8,:).';

% matlab imprementation
[PC, P] = partialcorr(X.');
pc1 = PC(1,2);
pc2 = partialcorr(x,y,z);
% regression version (original definition)
[b1,bint1,r1] = regress(x,[z, ones(sigLen,1)]);
[b2,bint2,r2] = regress(y,[z, ones(sigLen,1)]);
pc3 = (r1.'*r2) / (sqrt(r1.'*r1)*sqrt(r2.'*r2));
pc4 = corr(r1,r2);
% invert of cov matrix version
sigma = cov(X.');
B = inv(sigma);
pc5 = -B(1,2) / sqrt(B(1,1)*B(2,2));

PC3 = calcPartialCorrelation__(X);
Z = PC - PC3;
figure; clims = [-1 1]; imagesc(Z,clims); title('PC - PC3');

% synchronize signal 6 == 2, 7 invert 3
% this is rank down operation, then invert of cov matrix does not work!!
X(2,:) = X(6,:);
X(3,:) = 1 - X(7,:);

PC = calcPartialCorrelation(X); % calc PC - many NaN values

PC2 = calcPartialCorrelation_(X,0); % calc PC - this does not work well.

PC3 = calcPartialCorrelation__(X); % calc PC

PC4 = calcPLSPartialCorrelation(X); % calc PLS PC

% plot matrix
figure;
clims = [-1 1];
imagesc(PC,clims);
title('Partial Correlation (matlab implementation)');
colorbar;

% plot matrix
figure;
clims = [-1 1];
imagesc(PC2,clims);
title('Partial Correlation (invert of cov matrix)');
colorbar;

% plot matrix
figure;
clims = [-1 1];
imagesc(PC3,clims);
title('Partial Correlation (regression base)');
colorbar;

% plot matrix
figure;
clims = [-1 1];
imagesc(PC4,clims);
title('PLS Partial Correlation');
colorbar;

