% load experimental signals
subjectNum = 3;
roiNum = 52;
start = 4;

for idx=1:subjectNum
    expfile = ['data/ROI_Subject00' num2str(idx) '_Session001.mat'];

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

    outfName = ['data/marmoset-aneth-sample' num2str(idx) '-roi' num2str(roiNum) '.mat'];
    save(outfName, 'si', 'names');

    % check translated bold signal histgram
    figure; histogram(si(:));
    title([num2str(roiNum) 'ROI signal original histogram (' num2str(idx) ')']);

    % gaussian distribution to uniform like
    [Y, sig, m, maxsi, minsi] = gaussian2uniformSignal(si);
    figure; histogram(Y(:));
    title([num2str(roiNum) 'ROI signal converted histogram (' num2str(idx) ')']);
    
    % uniform to gaussian distribution (invert)
    si2 = uniform2gaussianSignal(Y, sig, m, maxsi, minsi);
    figure; histogram(si2(:));
    title([num2str(roiNum) 'ROI signal re-converted histogram (' num2str(idx) ')']);

    % diff between original and re-converted
    d = si - si2;
    err = sum(sum(d.*d))
end
