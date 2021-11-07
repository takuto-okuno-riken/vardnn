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

%% test pattern 1 -- just random
% regression version
PC3 = calcPartialCorrelation__(X);
Z = PC - PC3;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC3 : sum err=' num2str(nansum(abs(Z),'all'))]);

PC4 = calcPLSPartialCorrelation(X); % calc PLS PC
Z = PC - PC4;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC4 : sum err=' num2str(nansum(abs(Z),'all'))]);

PC5 = calcLassoPartialCorrelation(X, [], [], [], 0.01, 1); % calc Lasso PC
Z = PC - PC5;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC5 : sum err=' num2str(nansum(abs(Z),'all'))]);

[lambda, alpha, errMat] = estimateLassoParamsForPC(X, [], [], [], 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
PC5b = calcLassoPartialCorrelation(X, [], [], [], lambda, alpha); % calc Lasso PC
Z = PC - PC5b; nansum(abs(Z),'all')

PC6 = calcPcPartialCorrelation(X); % calc PCA+PC
Z = PC - PC6;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC6 : sum err=' num2str(nansum(abs(Z),'all'))]);

PC7 = calcSvPartialCorrelation(X); % calc SVR(linear)+PC
Z = PC - PC7;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC7 : sum err=' num2str(nansum(abs(Z),'all'))]);

PC8 = calcSvPartialCorrelation(X,[],[],[],'gaussian'); % calc SVR(gaussian)+PC
Z = PC - PC8;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC8 : sum err=' num2str(nansum(abs(Z),'all'))]);

PC9 = calcSvPartialCorrelation(X,[],[],[],'rbf'); % calc SVR(rbf)+PC
Z = PC - PC9;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC9 : sum err=' num2str(nansum(abs(Z),'all'))]);

%% test pattern 2 -- one step
exNum = 2;
exSignal = si(nodeNum+1:nodeNum+exNum,1:sigLen);
% control is all positive input
exControl = ones(nodeNum,exNum);
X = si(1:8, 1:sigLen);
X(2,2:end) = X(6,1:sigLen-1);
X(4,2:end) = X(6,1:sigLen-1); % caution! node 2 & 4 is Multicollinearity case (correlated)
X(5,2:end) = exSignal(1,1:sigLen-1);
%X(7,2:end) = exSignal(1,1:sigLen-1);

figure; PC = plotPartialCorrelation(X, exSignal, [], exControl, 1); % calc PC - many NaN values
PC3 = calcPartialCorrelation__(X, exSignal, [], exControl, 1);
figure; PC4 = plotPLSPartialCorrelation(X, exSignal, [], exControl, 1); % calc PLS PC
figure; PC5 = plotLassoPartialCorrelation(X, exSignal, [], exControl, 0.01, 1, 1); % calc Lasso PC
figure; PC6 = plotPcPartialCorrelation(X, exSignal, [], exControl, 1, 1); % calc PCA+PC
figure; PC7 = plotSvPartialCorrelation(X, exSignal, [], exControl, 'linear', 'auto', 1); % calc SVR(linber)+PC
figure; PC8 = plotSvPartialCorrelation(X, exSignal, [], exControl, 'gaussian', 'auto', 1); % calc SVR(gaussian)+PC
figure; PC9 = plotSvPartialCorrelation(X, exSignal, [], exControl, 'rbf', 'auto', 1); % calc SVR(rbf)+PC

Z = PC - PC3;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC3 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC4;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC4 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC5;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC5 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC6;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC6 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC7;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC7 : sum err=' num2str(nansum(abs(Z),'all'))]);
Z = PC - PC8;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC8 : sum err=' num2str(nansum(abs(Z),'all'))]);
Z = PC - PC9;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC9 : sum err=' num2str(nansum(abs(Z),'all'))]);


%% test pattern 3 -- one step
exNum = 8;
exSignal = si(nodeNum+1:nodeNum+exNum,1:sigLen);
% control is all positive input
exControl = eye(nodeNum,exNum);
X = si(1:8, 1:sigLen);
X(2,2:end) = X(6,1:sigLen-1);
X(4,2:end) = X(6,1:sigLen-1); % caution! node 2 & 4 is Multicollinearity case (correlated)
X(5,2:end) = exSignal(5,1:sigLen-1);
X(7,2:end) = exSignal(7,1:sigLen-1);

figure; PC = plotPartialCorrelation(X, exSignal, [], exControl, 1); % calc PC - many NaN values
PC3 = calcPartialCorrelation__(X, exSignal, [], exControl, 1);
figure; PC4 = plotPLSPartialCorrelation(X, exSignal, [], exControl, 1); % calc PLS PC
figure; PC5 = plotLassoPartialCorrelation(X, exSignal, [], exControl, 0.01, 1, 1); % calc Lasso PC
figure; PC6 = plotPcPartialCorrelation(X, exSignal, [], exControl, 1, 1); % calc PCA+PC
figure; PC7 = plotSvPartialCorrelation(X, exSignal, [], exControl, 'linear', 'auto', 1); % calc SVR(linber)+PC
figure; PC8 = plotSvPartialCorrelation(X, exSignal, [], exControl, 'gaussian', 'auto', 1); % calc SVR(gaussian)+PC
figure; PC9 = plotSvPartialCorrelation(X, exSignal, [], exControl, 'rbf', 'auto', 1); % calc SVR(rbf)+PC

Z = PC - PC3;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC3 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC4;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC4 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC5;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC5 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC6;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC6 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC - PC7;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC7 : sum err=' num2str(nansum(abs(Z),'all'))]);
Z = PC - PC8;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC8 : sum err=' num2str(nansum(abs(Z),'all'))]);
Z = PC - PC9;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC - PC9 : sum err=' num2str(nansum(abs(Z),'all'))]);

%% test pattern 4 -- rank down
% synchronize signal 6 == 2, 7 invert 3
% this is rank down operation, then invert of cov matrix does not work!!
X = si(1:8, 1:sigLen);
X(2,:) = X(6,:);
X(3,:) = 1 - X(7,:);

figure; PC = plotPartialCorrelation(X); % calc PC - many NaN values
PC2 = calcPartialCorrelation_(X,0); % calc PC - this does not work well.
PC3 = calcPartialCorrelation__(X); % calc PC
figure; PC4 = plotPLSPartialCorrelation(X); % calc PLS PC
figure; PC5 = plotLassoPartialCorrelation(X, [], [], [], 0.01, 1); % calc Lasso PC
figure; PC6 = plotPcPartialCorrelation(X, [], [], [], 1); % calc PCA+PC
figure; PC7 = plotSvPartialCorrelation(X, [], [], [], 'linear', 'auto', 1); % calc SVR(linber)+PC
figure; PC8 = plotSvPartialCorrelation(X, [], [], [], 'gaussian', 'auto', 1); % calc SVR(gaussian)+PC
figure; PC9 = plotSvPartialCorrelation(X, [], [], [], 'rbf', 'auto', 1); % calc SVR(rbf)+PC

[lambda, alpha, errMat] = estimateLassoParamsForPC(X, [], [], [], 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
PC5b = calcLassoPartialCorrelation(X, [], [], [], lambda, alpha); % calc Lasso PC
Z = PC3 - PC5b; nansum(abs(Z),'all')

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
Z = PC3 - PC4;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC3 - PC4 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC3 - PC5;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC3 - PC5 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC3 - PC6;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC3 - PC6 : sum err=' num2str(nansum(abs(Z),'all'))]);

Z = PC3 - PC7;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC3 - PC7 : sum err=' num2str(nansum(abs(Z),'all'))]);
Z = PC3 - PC8;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC3 - PC8 : sum err=' num2str(nansum(abs(Z),'all'))]);
Z = PC3 - PC9;
figure; clims = [-1 1]; imagesc(Z,clims); title(['PC3 - PC9 : sum err=' num2str(nansum(abs(Z),'all'))]);
