% load experimental signals
ROINUM = 132;
subjectNums = [32, 41, 8];
groups = {'ad', 'cn', 'mci'};

for k=1:length(groups)
    group = groups{k};
    subjectNum = subjectNums(k);
    outfName = ['data/ad-signal-' group '-roi' num2str(ROINUM) '.mat'];
    load(outfName);
    for i=1:subjectNum
        dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
        load(dlcmName); %, 'netDLCM', 'si', 'inSignal', 'inControl', 'mat', 'sig', 'c', 'maxsi', 'minsi');
        siorg = convert2InvSigmoidSignal(si, sig, m, maxsi, minsi);
        diff = siorg - signals{i};
        diffall = nansum(diff(:));
        disp([group '-' num2str(i) ' diff=' num2str(diffall)]);
        signals{i} = siorg;
    end
%    save(outfName, 'signals', 'roiNames');
end
