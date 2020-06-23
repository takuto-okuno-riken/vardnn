% load experimental signals
subjectNum = 3;
roiNum = 52;
start = 4;
sig_a = 0.5;

for idx=1:subjectNum
    expfile = ['test/ROI_Subject00' num2str(idx) '_Session001.mat'];

    %load conn ROI signal file
    load(expfile);

    regionNum = roiNum * 2;
    seqLen = size(data{1,start},1);
    names2 = names;

    si = zeros(regionNum,seqLen);
    names = cell(1,regionNum);

    for i=1:regionNum
        si(i,:) = data{1,start+(i-1)}.';
        names{1,i} = names2{1,start+(i-1)};
    end

    % save signal
    figure;
    plot(si.');
    title([num2str(roiNum) 'ROI signals (' num2str(idx) ')']);

    outfName = ['test/marmoset-aneth-sample' num2str(idx) '-roi' num2str(roiNum) '.mat'];
    save(outfName, 'si', 'names', 'sig_a');

    % check translated bold signal histgram
    si = bold2dnnSignal(si,sig_a);
    figure;
    histogram(si(:));
    title([num2str(roiNum) 'ROI signal histogram (' num2str(idx) ')']);
end
