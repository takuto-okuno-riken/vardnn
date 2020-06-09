nodeNum = 8;  % node number
sigLen = 100; % signal length

% generate random signals
X = rand(nodeNum, sigLen); 
Y = X; % keep original

% synchronize signal 6 == 2, 7 invert 3
X(2,:) = X(6,:);
X(3,:) = 1 - X(7,:);

FC = calcTimeShiftedCorrelation(X,0); % calc FC

% plot matrix
figure;
clims = [-1 1];
imagesc(FC,clims);
title('Sliding Functional Connectivity 0 : sync 6==2, 7 inv 3');
colorbar;

% synchronize signal 6 == 2, 7 invert 3
FC = calcTimeShiftedCorrelation(X,1); % calc FC

% plot matrix
figure;
clims = [-1 1];
imagesc(FC,clims);
title('Sliding Functional Connectivity 1 : sync 6==2, 7 inv 3');
colorbar;

% signal 6->2, 7-|3
X = Y;
X(2,2:sigLen) = X(6,1:sigLen-1);
X(3,2:sigLen) = 1 - X(7,1:sigLen-1);
FC = calcTimeShiftedCorrelation(X,1); % calc FC

% plot matrix
figure;
clims = [-1 1];
imagesc(FC,clims);
title('Sliding Functional Connectivity 1 : sync 6->2, 7-|3');
colorbar;

% signal 6->2, 7-|3
FC = calcTimeShiftedCorrelation(X,2); % calc FC

% plot matrix
figure;
clims = [-1 1];
imagesc(FC,clims);
title('Sliding Functional Connectivity 2 : sync 6->2, 7-|3');
colorbar;

% signal 6->2, 7-|3
FC = calcTimeShiftedCorrelation(X,3); % calc FC

% plot matrix
figure;
clims = [-1 1];
imagesc(FC,clims);
title('Sliding Functional Connectivity 3 : sync 6->2, 7-|3');
colorbar;
