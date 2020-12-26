% this function works after analyzeAlzheimerDLCM.m

function simulateAlzheimerDLCM
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');

    % simulate Cn signals
    [cnDLWs, smcnSignals] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw');
    [smcnDLs, ~, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlcm', 1);
    [smcnDLWs, meanSmcnDLW, stdSmcnDLW] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlw', 1);
    meanCnDLW = nanmean(cnDLWs,3);
    figure; cnsmcnDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcnDLW);

    % change Z score
%{
    cnDLWs = calcZScores(cnDLWs);
    adDLWs = calcZScores(adDLWs);
%}
    % plot correlation and cos similarity
    algNum = 1;
    meanCnDLW = nanmean(cnDLWs,3);
%    meanAdDLW = nanmean(adDLWs,3);
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    X = categorical({'cn-smcn'});
    figure; bar(X, cosSim);
    title('cos similarity between CN and SimCN by each algorithm');

    % normality test
%{
    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cnec', 'dlw');
%}
end

% ==================================================================================================================

function [ECs, simSignals] = simulateNodeSignals(signals, roiNames, group, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    outfName = ['results/adsim2-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    ECs = zeros(ROINUM, ROINUM, sbjNum);
    simSignals = cell(1, sbjNum);
    for i=1:sbjNum
        switch(algorithm)
        case 'dlw'
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            load(dlcmName);
            [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
            [Y, time] = simulateDlcmNetwork(si, inSignal, [], inControl, netDLCM);
            ec = calcDlcmEC(netDLCM, [], inControl);
        end
        ECs(:,:,i) = ec;
        simSignals{i} = Y;
    end
    save(outfName, 'ECs', 'simSignals', 'roiNames');
end


function [utestH, utestP, utestP2] = calculateAlzWilcoxonTest(control, target, roiNames, controlGroup, targetGroup, algorithm, force, testname, showfig)
    if nargin < 9
        showfig = 1;
    end
    if nargin < 8
        testname = 'ranksum';
    end
    if nargin < 7
        force = 0;
    end
    % constant value
    ROWNUM = size(control,1);
    COLNUM = size(control,2);

    outfName = ['results/adsim-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROWNUM) '-utest.mat'];
    if exist(outfName, 'file') && force == 0
        load(outfName);
    else
        utestH = nan(ROWNUM, COLNUM);
        utestP = nan(ROWNUM, COLNUM);
        utestP2 = nan(ROWNUM, COLNUM);
        for i=1:ROWNUM
            for j=1:COLNUM
                if isnan(control(i,j,1)), continue; end
                x = squeeze(control(i,j,:));
                y = squeeze(target(i,j,:));
                if length(x) > length(y)
                    x2 = x;
                    y2 = nan(length(x),1);
                    y2(1:length(y),1) = y;
                else
                    x2 = nan(length(y),1);
                    x2(1:length(x),1) = x;
                    y2 = y;
                end
                if isempty(find(~isnan(x2))) || isempty(find(~isnan(y2))), continue; end
                switch(testname)
                    case 'ranksum'
                        [p, h] = ranksum(x2,y2);
                    case 'ttest2'
                        [h, p] = ttest2(x2,y2);
                end
                %[p, h] = signrank(x2,y2);
                utestH(i,j) = h;
                utestP(i,j) = p;
                if h > 0 && nanmean(x) > nanmean(y)
                    utestP2(i,j) = p;
                end
            end
        end
        if force == 0
            save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames');
        end
    end
    % counting by source region and target region
    countSource = nansum(utestH,1);
    countTarget = nansum(utestH,2);
    save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames', 'countSource', 'countTarget');
    if showfig == 0, return; end
    
    load('test/colormap.mat')
    % U test result
    figure; 
    colormap(hvalmap);
    clims = [0,1];
    imagesc(utestH,clims);
    daspect([1 1 1]);
    title([controlGroup '-' targetGroup ' : ' algorithm ' : u test result']);
    colorbar;
    % U test p values
    utestP(isnan(utestP)) = 1;
    figure;
    colormap(pvalmap);
    clims = [0,1];
    imagesc(utestP, clims);
    daspect([1 1 1]);
    title([controlGroup '-' targetGroup ' : ' algorithm ' : u test p values']);
    colorbar;
end
