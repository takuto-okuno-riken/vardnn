nodeNum = 8;  % node number
sigLen = 200; % signal length
p = 3; % TE lag param

% generate random signals
%X = rand(nodeNum, sigLen); 

load('test/testTrain-rand500-uniform.mat');
X = si(1:8, 1:sigLen);

% set signal time lag 6->2, 6->4
X(2,3:end) = X(6,2:sigLen-1);
X(4,2:end) = X(6,1:sigLen-1);

%X(2,2:end) = X(6,1:sigLen-1);
%X(4,3:end) = X(2,2:sigLen-1);

TE = calcLinueTE(X, p); % calc transfer entropy index of lag |p|

% plot matrix
figure;
clims = [0 10];
imagesc(TE,clims);
title('Transfer Entropy (LINUE)');
colorbar;

