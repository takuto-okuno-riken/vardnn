nodeNum = 8;  % node number
sigLen = 100; % signal length

% generate random signals
X = rand(nodeNum, sigLen); 

% synchronize signal 6 == 2, 7 invert 3
X(2,:) = X(6,:);
X(3,:) = 1 - X(7,:);

FC = calcFunctionalConnectivity(X); % calc FC

% plot matrix
figure;
clims = [-1 1];
imagesc(FC,clims);
title('Functional Connectivity');
colorbar;

