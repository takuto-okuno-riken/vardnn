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

% find best lag number
maxLag = 10;
TE3 = zeros(nodeNum,nodeNum,maxLag);
nodeAIC = zeros(nodeNum,maxLag);
nodeBIC = zeros(nodeNum,maxLag);
for k=1:maxLag
    [TE3(:,:,k), h, P, ~, ~, ~, ~, nodeAIC(:,k), nodeBIC(:,k)] = calcLinueTE(X, k);
    h2 = double(P<0.05);
    h(isnan(h)) = 0;
    if isequal(h,h2)
        disp("h == P<0.05 !");
    else
        disp("error : h == P<0.05 !");
        return;
    end
end

% plot AIC and BIC of each node by each lag
figure; plot(nodeAIC.');
title('AIC of each node by each lag');
figure; plot(nodeBIC.');
title('BIC of each node by each lag');
[M,I]=min(nodeAIC,[],2);
lag = mode(I);
figure; histogram(I);
title('histogram of minimum AIC (lag) of each node');

% plot matrix
figure;
clims = [0 10];
imagesc(TE3(:,:,lag),clims); % actually, this result is not so good...
title('Transfer Entropy (LINUE) Index2');
colorbar;
