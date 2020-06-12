% set global random stream and shuffle it
myStream=RandStream('mt19937ar');
RandStream.setGlobalStream(myStream);
rng('shuffle');

% generate random signals
regionNum = 500;
seqLen = 2000;
si = rand(regionNum,seqLen);
plot(si.');

% generate random teaching signals
z = rand(1,seqLen);
figure;
plot(z);

save(['test/testTrain-rand' num2str(regionNum) '-uniform.mat'], 'si', 'z');
