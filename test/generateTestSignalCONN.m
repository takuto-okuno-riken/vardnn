% load experimental signals
idx = 3;
expfile = ['test/ROI_Subject00' num2str(idx) '_Session001.mat'];
load(expfile);

start = 4;
regionNum = 225 * 2;
seqLen = size(data{1,start},1);
names2 = names;

si = zeros(regionNum,seqLen);
names = cell(1,regionNum);

for i=1:regionNum
    si(i,:) = data{1,start+(i-1)}.';
    names{1,i} = names2{1,start+(i-1)};
end

figure;
plot(si.');

save(['test/marmoset-aneth-sample' num2str(idx) '-roi225.mat'], 'si', 'names');
