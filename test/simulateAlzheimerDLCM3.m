% this function works after analyzeAlzheimerDLCM2.m

function simulateAlzheimerDLCM3
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');

    cnSbjNum = length(cnSignals);
    adSbjNum = length(adSignals);
    nodeNum = size(cnSignals{1},1);
    maxLag = 5;

    % get DLCM-GC from original CN and AD (normal and recovery training)
    [cnDLs, meanCnDL, stdCnDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 1, 1, 1);
    [adDLs, meanAdDL, stdAdDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 1, 1, 1);
    sigCnDLs = (cnDLs - nanmean(cnDLs(:))) / nanstd(cnDLs(:),1);
    sigAdDLs = (adDLs - nanmean(adDLs(:))) / nanstd(adDLs(:),1);

    % simulate by algorithms
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [cnDLWs{j}, smcnSignals{j}, cnSubDLWs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn', 0, 0, j, 0);
        [adDLWs{j}, smadSignals{j}, adSubDLWs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad', 0, 0, j, 0);
        sigCnDLWs{j} = (cnDLWs{j} - nanmean(cnDLWs{j}(:))) / nanstd(cnDLWs{j}(:),1);
        sigAdDLWs{j} = (adDLWs{j} - nanmean(adDLWs{j}(:))) / nanstd(adDLWs{j}(:),1);
        meanCnDLW{j} = nanmean(cnDLWs{j},3);
        meanAdDLW{j} = nanmean(adDLWs{j},3);

        % DLCM(j) linear no exogenous 
        [cnDLW2s{j}, sm2cnSignals{j}, cnSubDLW2s{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn', 0, 0, j, 0, []);
        [adDLW2s{j}, sm2adSignals{j}, adSubDLW2s{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad', 0, 0, j, 0, []);
        sigCnDLW2s{j} = (cnDLW2s{j} - nanmean(cnDLW2s{j}(:))) / nanstd(cnDLW2s{j}(:),1);
        sigAdDLW2s{j} = (adDLW2s{j} - nanmean(adDLW2s{j}(:))) / nanstd(adDLW2s{j}(:),1);
        meanCnDLW2{j} = nanmean(cnDLW2s{j},3);
        meanAdDLW2{j} = nanmean(adDLW2s{j},3);

        % mvar(j) no exogenous 
        [cnMVARECs{j}, smmvcnSignals{j}, cnSubMVARECs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'mvarec', 'cn', 0, 0, j, 0);
        [adMVARECs{j}, smmvadSignals{j}, adSubMVARECs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'mvarec', 'ad', 0, 0, j, 0);
        sigCnMVARECs{j} = (cnMVARECs{j} - nanmean(cnMVARECs{j}(:))) / nanstd(cnMVARECs{j}(:),1);
        sigAdMVARECs{j} = (adMVARECs{j} - nanmean(adMVARECs{j}(:))) / nanstd(adMVARECs{j}(:),1);
        meanCnMVAREC{j} = nanmean(cnMVARECs{j},3);
        meanAdMVAREC{j} = nanmean(adMVARECs{j},3);
    end

    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [cnDLWs{j}, smcnSignals{j}, cnSubDLWs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn', 0, 0, j, 1);
        [adDLWs{j}, smadSignals{j}, adSubDLWs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad', 0, 0, j, 1);
        sigCnDLWs{j} = (cnDLWs{j} - nanmean(cnDLWs{j}(:))) / nanstd(cnDLWs{j}(:),1);
        sigAdDLWs{j} = (adDLWs{j} - nanmean(adDLWs{j}(:))) / nanstd(adDLWs{j}(:),1);
        meanCnDLW{j} = nanmean(cnDLWs{j},3);
        meanAdDLW{j} = nanmean(adDLWs{j},3);

        % DLCM(j) linear auto exogenous 
        [cnDLW2s{j}, sm2cnSignals{j}, cnSubDLW2s{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn', 0, 0, j, 1, []);
        [adDLW2s{j}, sm2adSignals{j}, adSubDLW2s{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad', 0, 0, j, 1, []);
        sigCnDLW2s{j} = (cnDLW2s{j} - nanmean(cnDLW2s{j}(:))) / nanstd(cnDLW2s{j}(:),1);
        sigAdDLW2s{j} = (adDLW2s{j} - nanmean(adDLW2s{j}(:))) / nanstd(adDLW2s{j}(:),1);
        meanCnDLW2{j} = nanmean(cnDLW2s{j},3);
        meanAdDLW2{j} = nanmean(adDLW2s{j},3);

        % mvar(j) auto exogenous 
        [cnMVARECs{j}, smmvcnSignals{j}, cnSubMVARECs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'mvarec', 'cn', 0, 0, j, 1);
        [adMVARECs{j}, smmvadSignals{j}, adSubMVARECs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'mvarec', 'ad', 0, 0, j, 1);
        sigCnDLWs{j} = (cnMVARECs{j} - nanmean(cnMVARECs{j}(:))) / nanstd(cnMVARECs{j}(:),1);
        sigAdDLWs{j} = (adMVARECs{j} - nanmean(adMVARECs{j}(:))) / nanstd(adMVARECs{j}(:),1);
        meanCnDLW{j} = nanmean(cnMVARECs{j},3);
        meanAdDLW{j} = nanmean(adMVARECs{j},3);
    end

    % --------------------------------------------------------------------------------------------------------------
    % check DLCM-EC and DLCM-GC of simulated CN and AD
    [smcnDLs, meanSmcnDL, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlcm', 1, 1, 1);
    [smcnDLWs, meanSmcnDLW, ~, smcnSubDLWs] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlw', 1, 1, 1);
    [smadDLs, meanSmadDL, ~] = calculateConnectivity(smadSignals, roiNames, 'smad', 'dlcm', 1, 1, 1);
    [smadDLWs, meanSmadDLW, ~, smadSubDLWs] = calculateConnectivity(smadSignals, roiNames, 'smad', 'dlw', 1, 1, 1);
    sigSmcnDLWs = (smcnDLWs - nanmean(smcnDLWs(:))) / nanstd(smcnDLWs(:),1);
    sigSmadDLWs = (smadDLWs - nanmean(smadDLWs(:))) / nanstd(smadDLWs(:),1);

    % check FC of simulated CN and AD
    [cnFCs, meanCnFC, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc', 1);
    [adFCs, meanAdFC, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc', 1);
    [smcnFCs, meanSmcnFC, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'fc', 1);
    [smadFCs, meanSmadFC, ~] = calculateConnectivity(smadSignals, roiNames, 'smad', 'fc', 1);
%    [smrccnFCs, meanSmrccnFC, ~] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'fc', 1);
%    [smrcadFCs, meanSmrcadFC, ~] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'fc', 1);
    [smcn2FCs, meanSmcn2FC, ~] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'fc', 1);
    [smad2FCs, meanSmad2FC, ~] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'fc', 1);
    [smcn6FCs, meanSmcn6FC, ~] = calculateConnectivity(smcn6Signals, roiNames, 'smcn6', 'fc', 1);
    [smcn7FCs, meanSmcn7FC, ~] = calculateConnectivity(smcn7Signals, roiNames, 'smcn7', 'fc', 1);
    [smcn8FCs, meanSmcn8FC, ~] = calculateConnectivity(smcn8Signals, roiNames, 'smcn8', 'fc', 1);
    [smcn9FCs, meanSmcn9FC, ~] = calculateConnectivity(smcn9Signals, roiNames, 'smcn9', 'fc', 1);
    [smad7FCs, meanSmad7FC, ~] = calculateConnectivity(smad7Signals, roiNames, 'smad7', 'fc', 1);
    [smad8FCs, meanSmad8FC, ~] = calculateConnectivity(smad8Signals, roiNames, 'smad8', 'fc', 1);
    [wtcnFCs, meanWtcnFC, ~] = calculateConnectivity(wtcnSignals, roiNames, 'wtcn', 'fc', 1);
%    [wtcn2FCs, meanWtcn2FC, ~] = calculateConnectivity(wtcn2Signals, roiNames, 'wtcn2', 'fc', 1);
    [wtcn3FCs, meanWtcn3FC, ~] = calculateConnectivity(wtcn3Signals, roiNames, 'wtcn3', 'fc', 1);
    [wtcn4FCs, meanWtcn4FC, ~] = calculateConnectivity(wtcn4Signals, roiNames, 'wtcn4', 'fc', 1);
    smmvcnFCs     = cell(lagMax,2);
    meanSmmvcnFC  = cell(lagMax,2);
    smmvadFCs     = cell(lagMax,2);
    meanSmmvadFC  = cell(lagMax,2);
    for i=1:lagMax
        [smmvcnFCs{i,1},  meanSmmvcnFC{i,1}, ~] = calculateConnectivity(smmvcnSignals{i,1}, roiNames, ['smmvcn' num2str(i)], 'fc', 1);
        [smmvadFCs{i,1},  meanSmmvadFC{i,1}, ~] = calculateConnectivity(smmvadSignals{i,1}, roiNames, ['smmvad' num2str(i)], 'fc', 1);        
        [smmvcnFCs{i,2},  meanSmmvcnFC{i,2}, ~] = calculateConnectivity(smmvcnSignals{i,2}, roiNames, ['smmvcn2' num2str(i)], 'fc', 1);
        [smmvadFCs{i,2},  meanSmmvadFC{i,2}, ~] = calculateConnectivity(smmvadSignals{i,2}, roiNames, ['smmvad2' num2str(i)], 'fc', 1);
    end

%{
    figure; cnsmcnFCr = plotTwoSignalsCorrelation(meanCnFC, meanSmcnFC);
    figure; adsmadFCr = plotTwoSignalsCorrelation(meanAdFC, meanSmadFC);
    figure; cnsmrccnFCr = plotTwoSignalsCorrelation(meanCnFC, meanSmrccnFC);
    figure; adsmrcadFCr = plotTwoSignalsCorrelation(meanAdFC, meanSmrcadFC);
    figure; cnsmcn2FCr = plotTwoSignalsCorrelation(meanCnFC, meanSmcn2FC);
    figure; adsmad2FCr = plotTwoSignalsCorrelation(meanAdFC, meanSmad2FC);
%}
    
    % --------------------------------------------------------------------------------------------------------------
    % plot correlation and cos similarity
    algNum = 30;
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    cosSim(2) = getCosSimilarity(meanCnDLW, meanSmcn2DLW);
    cosSim(3) = getCosSimilarity(meanCnDLW, meanSmcn6DLW);
    cosSim(4) = getCosSimilarity(meanCnDLW, meanSmcn7DLW);
    cosSim(5) = getCosSimilarity(meanCnDLW, meanSmcn8DLW);
    cosSim(6) = getCosSimilarity(meanCnDLW, meanSmcn9DLW);
    cosSim(7) = getCosSimilarity(meanCnDLW, meanWtcnDLW);
%    cosSim(8) = getCosSimilarity(meanCnDLW, meanWtcn2DLW);
    cosSim(8) = getCosSimilarity(meanCnDLW, meanWtcn3DLW);
    cosSim(9) = getCosSimilarity(meanCnDLW, meanWtcn4DLW);
    cosSim(11) = getCosSimilarity(meanCnDL, meanSmcnDL);
    cosSim(12) = getCosSimilarity(meanCnDL, meanSmcn2DL);
    cosSim(13) = getCosSimilarity(meanCnDL, meanSmcn7DL);
    cosSim(14) = getCosSimilarity(meanCnDL, meanSmcn8DL);
    cosSim(15) = getCosSimilarity(meanCnDL, meanSmcn9DL);
    cosSim(16) = getCosSimilarity(meanCnDL, meanWtcnDL);
%    cosSim(17) = getCosSimilarity(meanCnDL, meanWtcn2DL);
    cosSim(17) = getCosSimilarity(meanCnDL, meanWtcn3DL);
    cosSim(18) = getCosSimilarity(meanCnDL, meanWtcn4DL);
    cosSim(21) = getCosSimilarity(meanCnFC, meanSmcnFC);
    cosSim(22) = getCosSimilarity(meanCnFC, meanSmcn2FC);
    cosSim(23) = getCosSimilarity(meanCnFC, meanSmcn6FC);
    cosSim(24) = getCosSimilarity(meanCnFC, meanSmcn7FC);
    cosSim(25) = getCosSimilarity(meanCnFC, meanSmcn8FC);
    cosSim(26) = getCosSimilarity(meanCnFC, meanSmcn9FC);
    cosSim(27) = getCosSimilarity(meanCnFC, meanWtcnFC);
%    cosSim(28) = getCosSimilarity(meanCnFC, meanWtcn2FC);
    cosSim(28) = getCosSimilarity(meanCnFC, meanWtcn3FC);
    cosSim(29) = getCosSimilarity(meanCnFC, meanWtcn4FC);
%    X = categorical({'dlec-cn-smcn','dlec-cn-smcn2','dlec-cn-smcn6','dlec-cn-smcn7',...
%        'dlgc-cn-smcn','dlgc-cn-smcn2','dlgc-cn-smcn7',...
%        'fc-cn-smcn','fc-cn-smcn2','fc-cn-smcn6','fc-cn-smcn7'});
    figure; bar(cosSim); title('cos similarity between mean CN matrix and SimCN by each algorithm');

    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanAdDLW, meanSmadDLW);
    cosSim(2) = getCosSimilarity(meanAdDLW, meanSmad2DLW);
    cosSim(3) = getCosSimilarity(meanAdDLW, meanSmad7DLW);
    cosSim(4) = getCosSimilarity(meanAdDLW, meanSmad8DLW);
    cosSim(11) = getCosSimilarity(meanAdDL, meanSmadDL);
    cosSim(12) = getCosSimilarity(meanAdDL, meanSmad2DL);
    cosSim(13) = getCosSimilarity(meanAdDL, meanSmad7DL);
    cosSim(14) = getCosSimilarity(meanAdDL, meanSmad8DL);
    cosSim(21) = getCosSimilarity(meanAdFC, meanSmadFC);
    cosSim(22) = getCosSimilarity(meanAdFC, meanSmad2FC);
    cosSim(23) = getCosSimilarity(meanAdFC, meanSmad7FC);
    cosSim(24) = getCosSimilarity(meanAdFC, meanSmad8FC);
%    X = categorical({'dlec-ad-smad','dlec-ad-smad2','dlgc-ad-smad','dlgc-ad-smad2','fc-ad-smad','fc-ad-smad2','',''});
    figure; bar(cosSim); title('cos similarity between mean AD matrix and SimAD by each algorithm');
%{
    % plot correlation and cos similarity (all)
    algNum = 24;
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    cosSim(2) = getCosSimilarity(meanAdDLW, meanSmadDLW);
    cosSim(3) = getCosSimilarity(meanAdDLW, meanSmcnDLW);
    cosSim(4) = getCosSimilarity(meanCnDLW, meanSmadDLW);
    cosSim(5) = getCosSimilarity(meanCnDL, meanSmcnDL);
    cosSim(6) = getCosSimilarity(meanAdDL, meanSmadDL);
    cosSim(7) = getCosSimilarity(meanAdDL, meanSmcnDL);
    cosSim(8) = getCosSimilarity(meanCnDL, meanSmadDL);
    cosSim(9) = getCosSimilarity(meanCnFC, meanSmcnFC);
    cosSim(10) = getCosSimilarity(meanAdFC, meanSmadFC);
    cosSim(11) = getCosSimilarity(meanAdFC, meanSmcnFC);
    cosSim(12) = getCosSimilarity(meanCnFC, meanSmadFC);
%    cosSim(13) = getCosSimilarity(meanCnDLW, meanSmrccnDLW);
%    cosSim(14) = getCosSimilarity(meanAdDLW, meanSmrcadDLW);
%    cosSim(15) = getCosSimilarity(meanAdDLW, meanSmrccnDLW);
%    cosSim(16) = getCosSimilarity(meanCnDLW, meanSmrcadDLW);
%    cosSim(17) = getCosSimilarity(meanCnFC, meanSmrccnFC);
%    cosSim(18) = getCosSimilarity(meanAdFC, meanSmrcadFC);
%    cosSim(19) = getCosSimilarity(meanAdFC, meanSmrccnFC);
%    cosSim(20) = getCosSimilarity(meanCnFC, meanSmrcadFC);
    cosSim(21) = getCosSimilarity(meanCnDLW, meanSmcn2DLW);
    cosSim(22) = getCosSimilarity(meanAdDLW, meanSmad2DLW);
    cosSim(23) = getCosSimilarity(meanAdDLW, meanSmcn2DLW);
    cosSim(24) = getCosSimilarity(meanCnDLW, meanSmad2DLW);
    X = categorical({'cn-smcn','ad-smad','ad-smcn','cn-smad','cn-smcn-dlgc','ad-smad-dlgc','ad-smcn-dlgc','cn-smad-dlgc',...
        'cn-smcn-fc','ad-smad-fc','ad-smcn-fc','cn-smad-fc','cn-smrccn','ad-smrcad','ad-smrccn','cn-smrcad',...
        'cn-smrccn-fc','ad-smrcad-fc','ad-smrccn-fc','cn-smrcad-fc','cn-smcn2','ad-smad2','ad-smcn2','cn-smad2',});
    figure; bar(X, cosSim);
    title('cos similarity between CN and SimCN by each algorithm');
%}
    % plot box-and-whisker plot of cos similarity between mean ec matrix and each subject ec
    cosSims = nan(cnSbjNum,30);
    for i=1:cnSbjNum
        cosSims(i,1) = getCosSimilarity(meanCnDLW, cnDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(meanCnDLW, smcnDLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(meanCnDLW, smcn2DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(meanCnDLW, smcn6DLWs(:,:,i));
        cosSims(i,5) = getCosSimilarity(meanCnDLW, smcn7DLWs(:,:,i));
        cosSims(i,6) = getCosSimilarity(meanCnDLW, smcn8DLWs(:,:,i));
        cosSims(i,7) = getCosSimilarity(meanCnDLW, smcn9DLWs(:,:,i));
        cosSims(i,8) = getCosSimilarity(meanCnDLW, wtcnDLWs(:,:,i));
%        cosSims(i,9) = getCosSimilarity(meanCnDLW, wtcn2DLWs(:,:,i));
        cosSims(i,9) = getCosSimilarity(meanCnDLW, wtcn3DLWs(:,:,i));
        cosSims(i,10) = getCosSimilarity(meanCnDLW, wtcn4DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(meanCnDL, cnDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(meanCnDL, smcnDLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(meanCnDL, smcn2DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(meanCnDL, smcn7DLs(:,:,i));
        cosSims(i,15) = getCosSimilarity(meanCnDL, smcn8DLs(:,:,i));
        cosSims(i,16) = getCosSimilarity(meanCnDL, smcn9DLs(:,:,i));
        cosSims(i,17) = getCosSimilarity(meanCnDL, wtcnDLs(:,:,i));
%        cosSims(i,18) = getCosSimilarity(meanCnDL, wtcn2DLs(:,:,i));
        cosSims(i,18) = getCosSimilarity(meanCnDL, wtcn3DLs(:,:,i));
        cosSims(i,19) = getCosSimilarity(meanCnDL, wtcn4DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(meanCnFC, cnFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(meanCnFC, smcnFCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(meanCnFC, smcn2FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(meanCnFC, smcn6FCs(:,:,i));
        cosSims(i,25) = getCosSimilarity(meanCnFC, smcn7FCs(:,:,i));
        cosSims(i,26) = getCosSimilarity(meanCnFC, smcn8FCs(:,:,i));
        cosSims(i,27) = getCosSimilarity(meanCnFC, smcn9FCs(:,:,i));
        cosSims(i,28) = getCosSimilarity(meanCnFC, wtcnFCs(:,:,i));
%        cosSims(i,29) = getCosSimilarity(meanCnFC, wtcn2FCs(:,:,i));
        cosSims(i,29) = getCosSimilarity(meanCnFC, wtcn3FCs(:,:,i));
        cosSims(i,30) = getCosSimilarity(meanCnFC, wtcn4FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:30));
    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
    [sortCosCnDLW,idxCosCnDLW] = sort(cosSims(:,1),'descend');
    [sortCosCnFC,idxCosCnFC] = sort(cosSims(:,21),'descend');

    cosSims = nan(adSbjNum,30);
    for i=1:adSbjNum
        cosSims(i,1) = getCosSimilarity(meanAdDLW, adDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(meanAdDLW, smadDLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(meanAdDLW, smad2DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(meanAdDLW, smad7DLWs(:,:,i));
        cosSims(i,5) = getCosSimilarity(meanAdDLW, smad8DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(meanAdDL, adDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(meanAdDL, smadDLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(meanAdDL, smad2DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(meanAdDL, smad7DLs(:,:,i));
        cosSims(i,15) = getCosSimilarity(meanAdDL, smad8DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(meanAdFC, adFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(meanAdFC, smadFCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(meanAdFC, smad2FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(meanAdFC, smad7FCs(:,:,i));
        cosSims(i,25) = getCosSimilarity(meanAdFC, smad8FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:30));
    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
    [sortCosAdDLW,idxCosAdDLW] = sort(cosSims(:,1),'descend');
    [sortCosAdFC,idxCosAdFC] = sort(cosSims(:,21),'descend');

    % plot correlation between original EC matrix and simulated signals EC matrix
    for i=1:1 %cnSbjNum
        % calc & plot correlation of original Zi vs simulating Zi
        [corrZi, corrZij] = calcCorrelationZiZij(cnSubDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum);
        plotCorrelationZiZij(cnDLWs(:,:,i), cnSubDLWs(:,:,i), smcnDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum, ['sbj' num2str(i)], 'original', 'simulating');
    end

    % plot box-and-whisker plot of cos similarity between original ec matrix and simulated signals ec matrix
    cosSims = nan(cnSbjNum,30);
    for i=1:cnSbjNum
        cosSims(i,1) = getCosSimilarity(cnDLWs(:,:,i), smcnDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(cnDLWs(:,:,i), smcn2DLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(cnDLWs(:,:,i), smcn6DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(cnDLWs(:,:,i), smcn7DLWs(:,:,i));
        cosSims(i,5) = getCosSimilarity(cnDLWs(:,:,i), smcn8DLWs(:,:,i));
        cosSims(i,6) = getCosSimilarity(cnDLWs(:,:,i), smcn9DLWs(:,:,i));
        cosSims(i,7) = getCosSimilarity(cnDLWs(:,:,i), wtcnDLWs(:,:,i));
%        cosSims(i,8) = getCosSimilarity(cnDLWs(:,:,i), wtcn2DLWs(:,:,i));
        cosSims(i,8) = getCosSimilarity(cnDLWs(:,:,i), wtcn3DLWs(:,:,i));
        cosSims(i,9) = getCosSimilarity(cnDLWs(:,:,i), wtcn4DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(cnDLs(:,:,i), smcnDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(cnDLs(:,:,i), smcn2DLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(cnDLs(:,:,i), smcn7DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(cnDLs(:,:,i), smcn8DLs(:,:,i));
        cosSims(i,15) = getCosSimilarity(cnDLs(:,:,i), smcn9DLs(:,:,i));
        cosSims(i,16) = getCosSimilarity(cnDLs(:,:,i), wtcnDLs(:,:,i));
%        cosSims(i,17) = getCosSimilarity(cnDLs(:,:,i), wtcn2DLs(:,:,i));
        cosSims(i,17) = getCosSimilarity(cnDLs(:,:,i), wtcn3DLs(:,:,i));
        cosSims(i,18) = getCosSimilarity(cnDLs(:,:,i), wtcn4DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(cnFCs(:,:,i), smcnFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(cnFCs(:,:,i), smcn2FCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(cnFCs(:,:,i), smcn6FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(cnFCs(:,:,i), smcn7FCs(:,:,i));
        cosSims(i,25) = getCosSimilarity(cnFCs(:,:,i), smcn8FCs(:,:,i));
        cosSims(i,26) = getCosSimilarity(cnFCs(:,:,i), smcn9FCs(:,:,i));
        cosSims(i,27) = getCosSimilarity(cnFCs(:,:,i), wtcnFCs(:,:,i));
%        cosSims(i,28) = getCosSimilarity(cnFCs(:,:,i), wtcn2FCs(:,:,i));
        cosSims(i,28) = getCosSimilarity(cnFCs(:,:,i), wtcn3FCs(:,:,i));
        cosSims(i,29) = getCosSimilarity(cnFCs(:,:,i), wtcn4FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    cosSims = nan(cnSbjNum,60);
    for i=1:cnSbjNum
        for j=1:10
            cosSims(i,j) = getCosSimilarity(cnDLWs(:,:,i), smmvcnDLWs{j,1}(:,:,i));
            cosSims(i,10+j) = getCosSimilarity(cnDLWs(:,:,i), smmvcnDLWs{j,2}(:,:,i));
            cosSims(i,20+j) = getCosSimilarity(cnDLs(:,:,i), smmvcnDLs{j,1}(:,:,i));
            cosSims(i,30+j) = getCosSimilarity(cnDLs(:,:,i), smmvcnDLs{j,2}(:,:,i));
            cosSims(i,40+j) = getCosSimilarity(cnFCs(:,:,i), smmvcnFCs{j,1}(:,:,i));
            cosSims(i,50+j) = getCosSimilarity(cnFCs(:,:,i), smmvcnFCs{j,2}(:,:,i));
        end
    end
    figure; boxplot(cosSims);
    cosSims = nan(cnSbjNum,30);
    for i=1:adSbjNum
        cosSims(i,1) = getCosSimilarity(adDLWs(:,:,i), smadDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(adDLWs(:,:,i), smad2DLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(adDLWs(:,:,i), smad7DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(adDLWs(:,:,i), smad8DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(adDLs(:,:,i), smadDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(adDLs(:,:,i), smad2DLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(adDLs(:,:,i), smad7DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(adDLs(:,:,i), smad8DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(adFCs(:,:,i), smadFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(adFCs(:,:,i), smad2FCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(adFCs(:,:,i), smad7FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(adFCs(:,:,i), smad8FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    
    % change Z score
%{
    cnDLWs = calcZScores(cnDLWs);
    adDLWs = calcZScores(adDLWs);
%}

    % normality test
%{
    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cnec', 'dlw');
%}
    % compalizon test (Wilcoxon, Mann?Whitney U test)
    [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcnFCs, roiNames, 'cn', 'smcn', 'fc');
    [adsmadFCsUt, adsmadFCsUtP, adsmadFCsUtP2] = calculateAlzWilcoxonTest(adFCs, smadFCs, roiNames, 'ad', 'smad', 'fc');
%    [cnsmrccnFCsUt, cnsmrccnFCsUtP, cnsmrccnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smrccnFCs, roiNames, 'cn', 'smrccn', 'fc');
%    [adsmrcadFCsUt, adsmrcadFCsUtP, adsmrcadFCsUtP2] = calculateAlzWilcoxonTest(adFCs, smrcadFCs, roiNames, 'ad', 'smrcad', 'fc');
    [cnsmcn2FCsUt, cnsmcn2FCsUtP, cnsmcn2FCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcn2FCs, roiNames, 'cn', 'smcn2', 'fc');
    [adsmad2FCsUt, adsmad2FCsUtP, adsmad2FCsUtP2] = calculateAlzWilcoxonTest(adFCs, smad2FCs, roiNames, 'ad', 'smad2', 'fc');
    [cnsmcn6FCsUt, cnsmcn6FCsUtP, cnsmcn6FCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcn6FCs, roiNames, 'cn', 'smcn6', 'fc');
    [cnsmcn7FCsUt, cnsmcn7FCsUtP, cnsmcn7FCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcn7FCs, roiNames, 'cn', 'smcn7', 'fc');
%    [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, wtcnFCs, roiNames, 'cn', 'wtcn', 'fc');
    [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs, smcnDLs, roiNames, 'cn', 'smcn', 'dlcm');
    [adsmadDLsUt, adsmadDLsUtP, adsmadDLsUtP2] = calculateAlzWilcoxonTest(adDLs, smadDLs, roiNames, 'ad', 'smad', 'dlcm');
    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smcnDLWs, roiNames, 'cn', 'smcn', 'dlw');
    [adsmadDLWsUt, adsmadDLWsUtP, adsmadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smadDLWs, roiNames, 'ad', 'smad', 'dlw');
    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smcn7DLWs, roiNames, 'cn', 'smcn7', 'dlw');
%    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, wtcnDLWs, roiNames, 'cn', 'wtcn', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigCnDLWs, sigSmcnDLWs, roiNames, 'sigcn', 'sigsmcn', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigAdDLWs, sigSmadDLWs, roiNames, 'sigad', 'sigsmad', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigCnDLWs, sigSmcn7DLWs, roiNames, 'sigcn', 'sigsmcn7', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigAdDLWs, sigSmad7DLWs, roiNames, 'sigad', 'sigsmad7', 'dlw');
%    [cnsmrccnDLWsUt, cnsmrccnDLWsUtP, cnsmrccnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smrccnDLWs, roiNames, 'cn', 'smrccn', 'dlw');
%    [adsmrcadDLWsUt, adsmrcadDLWsUtP, adsmrcadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smrcadDLWs, roiNames, 'ad', 'smrcad', 'dlw');
end

% ==================================================================================================================

function plotCorrelationZiZij(EC, subEC, smEC, smSubEC, nodeNum, prefix, orig, sim)
    if ~isempty(subEC)
        figure; hold on; plot([0.5 1.2], [0.5 1.2],':','Color',[0.5 0.5 0.5]);
        for i=1:nodeNum
            X = subEC(i,2:end);
            Y = smSubEC(i,2:end); Y(i) = nan;
            plotTwoSignalsCorrelation(X, Y, [0.3+0.07*mod(i,10) 0.3+0.07*ceil(mod(i,100)/10) 0.6+0.2*ceil(i/100)]);
        end
        for i=1:nodeNum
            plotTwoSignalsCorrelation(subEC(i,1), smSubEC(i,1), [0.4+0.2*ceil(i/100) 0.1*mod(i,5) 0.1*ceil(mod(i,50)/10)], 'd', 8);
        end
        hold off; title([prefix ' Zij corr: ' orig ' vs ' sim]);
    end
    if ~isempty(EC)
        figure; hold on; plot([0 0.5], [0 0.5],':','Color',[0.5 0.5 0.5]);
        for i=1:nodeNum
            plotTwoSignalsCorrelation(EC(i,:), smEC(i,:), [0.4+0.2*ceil(i/100) 0.1*mod(i,5) 0.1*ceil(mod(i,50)/10)], 'x', 5);
        end
        hold off; title([prefix ' Zi corr: ' orig ' vs ' sim]);
    end
end

function [ZiCorr, ZijCorr] = calcCorrelationZiZij(subEC, smSubEC, nodeNum)        
    ZiCorr = corr2(subEC(:,1), smSubEC(:,1));
    ZijCorr = nan(nodeNum,1);
    for i=1:nodeNum
        X = subEC(i,2:end);
        Y = smSubEC(i,2:end); Y(i) = nan;
        ZijCorr(i) = corr2(X(~isnan(X)), Y(~isnan(Y)));
    end
end

function [ECs, simSignals, subECs] = simulateNodeSignals(signals, roiNames, group, algorithm, orgGroup, isRaw, inSiRange, lags, isAutoExo, activateFunc)
    if nargin < 10, activateFunc = @reluLayer; end
    if nargin < 9, isAutoExo = 0; end
    if nargin < 8, lags = 1; end
    if nargin < 7, inSiRange = 0; end
    if nargin < 6, isRaw = 0; end

    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    lagpat = ["gc","pgc","te","tsfc","tsfca","mvarec","dlcm","dlw"];
    if lags>1 && contains(algorithm,lagpat), lagStr=num2str(lags); else lagStr=''; end
    if isAutoExo>0, exoStr='_ex'; else exoStr=''; end
    if isempty(activateFunc), linStr='_lin'; else linStr=''; end
    outfName = ['results/adsim3-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 11;

    if NumProcessors > 1
        try
            disp('Destroing any existance matlab pool session');
            parpool('close');
        catch
            disp('No matlab pool session found');
        end
        parpool(NumProcessors);
    end
        
    ECs = zeros(ROINUM, ROINUM, sbjNum);
    subECs = zeros(ROINUM, ROINUM+1, sbjNum);
    simSignals = cell(1, sbjNum);
    parfor i=1:sbjNum    % for parallel processing
%    for i=1:sbjNum
        switch(algorithm)
        case 'dlw'
            dlcmName = ['results/ad-dlcm' lagStr exoStr linStr '-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            if isRaw
                si = signals{i};
                sig=0; c=0; maxsi=0; minsi=0;
            else
                [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
            end
            if inSiRange > 0
                exSignal = f.exSignal(:,inSiRange);
            else
                exSignal = f.exSignal;
            end
            [Y, time] = simulateDlcmNetwork(si, exSignal, [], f.exControl, f.netDLCM);
            [ec, subECs(:,:,i)] = calcDlcmEC(f.netDLCM, [], f.exControl);
        case 'mvarec'
            dlcmName = ['results/ad-dlcm' lagStr exoStr linStr '-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            if isRaw
                si = signals{i};
                sig=0; c=0; maxsi=0; minsi=0;
            else
                [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
            end
            if inSiRange > 0
                exSignal = f.exSignal(:,inSiRange);
            else
                exSignal = f.exSignal;
            end
            % simulate LAR network with 1st frame & exogenous input signal
            netMVAR = initMvarNetwork(si, exSignal, [], f.exControl, lags);
            [Y, time] = simulateMvarNetwork(si, exSignal, [], f.exControl, netMVAR);
            [ec, subECs(:,:,i)] = calcDlcmEC(f.netDLCM, [], f.exControl); % not used this
        end
        ECs(:,:,i) = ec;
        simSignals{i} = Y;
    end
    save(outfName, 'ECs', 'simSignals', 'roiNames', 'subECs');

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
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
