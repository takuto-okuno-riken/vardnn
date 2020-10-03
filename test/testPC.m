nodeNum = 8;  % node number
sigLen = 100; % signal length

% generate random signals
X = rand(nodeNum, sigLen); 

% synchronize signal 6 == 2, 7 invert 3
X(2,:) = X(6,:);
X(3,:) = 1 - X(7,:);

PC = calcPartialCorrelation(X); % calc PC

% plot matrix
figure;
clims = [-1 1];
imagesc(PC,clims);
title('Partial Correlation');
colorbar;

