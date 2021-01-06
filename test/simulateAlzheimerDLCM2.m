% this function works after analyzeAlzheimerDLCM.m

function simulateAlzheimerDLCM2
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
    
    % get DLCM-GC from original CN and AD (normal and recovery training)
    [cnDLs, meanCnDL, stdCnDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 1);
    [adDLs, meanAdDL, stdAdDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 1);
    sigCnDLs = (cnDLs - nanmean(cnDLs(:))) / nanstd(cnDLs(:),1);
    sigAdDLs = (adDLs - nanmean(adDLs(:))) / nanstd(adDLs(:),1);

%    [rccnDLs, meanRcCnDL, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcmrc', 1); % do recovery training
%    [rcadDLs, meanRcAdDL, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcmrc', 1); % do recovery training

    % simulate CN & AD signals from first frame
    [cnDLWs, smcnSignals, cnSubDLWs] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn');
    [adDLWs, smadSignals, adSubDLWs] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad');
    sigCnDLWs = (cnDLWs - nanmean(cnDLWs(:))) / nanstd(cnDLWs(:),1);
    sigAdDLWs = (adDLWs - nanmean(adDLWs(:))) / nanstd(adDLWs(:),1);
    meanCnDLW = nanmean(cnDLWs,3);
    meanAdDLW = nanmean(adDLWs,3);

    % simulate recovery trained CN & AD signals from first frame
%    [rccnDLWs, smrccnSignals] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlwrc', 'cn');
%    [rcadDLWs, smrcadSignals] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlwrc', 'ad');

    % simulate CN & AD signals from last-1 frame
    cn2Signals = cnSignals;
    ad2Signals = adSignals;
    for i=1:cnSbjNum
        cn2Signals{i}(:,1) = cn2Signals{i}(:,end-1);
        cn2Signals{i}(:,2) = cn2Signals{i}(:,end);
    end
    for i=1:adSbjNum
        ad2Signals{i}(:,1) = ad2Signals{i}(:,end-1);
        ad2Signals{i}(:,2) = ad2Signals{i}(:,end);
    end
    [~, smcn2Signals] = simulateNodeSignals(cn2Signals, roiNames, 'cn2', 'dlw', 'cn');
    [~, smad2Signals] = simulateNodeSignals(ad2Signals, roiNames, 'ad2', 'dlw', 'ad');

    % expanding amplitude of simulated CN & AD signals from last-1 frame (type1)
    % -- no effect
%{
    smcn3Signals = smcn2Signals;
    smad3Signals = smad2Signals;
    for i=1:cnSbjNum
        smcn3Signals{i} = expandAmplitude(smcn2Signals{i}, 1.3);
    end
    for i=1:adSbjNum
        smad3Signals{i} = expandAmplitude(smad2Signals{i}, 1.3);
    end
%}
    % expanding amplitude of simulated CN & AD signals from last-1 frame (type2)
    % finding best rate -- result showed bigger k is better (k=10 is best)
%{
    smcn4Signals = smcn2Signals;
    smad4Signals = smad2Signals;
    cnsmcn4DLWrs = zeros(1,10);
    for k=10:10
        for i=1:cnSbjNum
            smcn4Signals{i} = expandAmplitude2(smcn2Signals{i}, 1+0.5*k);
        end
        for i=1:adSbjNum
            smad4Signals{i} = expandAmplitude2(smad2Signals{i}, 1+0.5*k);
        end
        name = ['smcn4_' num2str(k)];
        [smcn4DLs, meanSmcn4DL, ~] = calculateConnectivity(smcn4Signals, roiNames, name, 'dlcm', 1);
        [smcn4DLWs, meanSmcn4DLW, ~, smcn4SubDLWs] = calculateConnectivity(smcn4Signals, roiNames, name, 'dlw', 1);
        figure; cnsmcn4DLWrs(k) = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn4DLW);
    end
%}
    % re-train CN & AD signals with expanded EC amplitude (type3)
    % -- no effect. trained DLCM network showed bad simulation signals
    S2 = ones(nodeNum, nodeNum+1);
    S2(:,2:end) = S2(:, 2:end) - eye(nodeNum);
    IS2 = ones(nodeNum, nodeNum+1);
    %[smcn5DLWs, smcn5Signals] = retrainDLCMAndECmultiPattern(cnSignals, cnDLWs, cnSubDLWs, S2, IS2, roiNames, 'cn5');

    % --------------------------------------------------------------------------------------------------------------
    % check DLCM-EC and DLCM-GC of simulated CN and AD
    [smcnDLs, meanSmcnDL, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlcm', 1);
    [smcnDLWs, meanSmcnDLW, ~, smcnSubDLWs] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlw', 1);
    [smadDLs, meanSmadDL, ~] = calculateConnectivity(smadSignals, roiNames, 'smad', 'dlcm', 1);
    [smadDLWs, meanSmadDLW, ~, smadSubDLWs] = calculateConnectivity(smadSignals, roiNames, 'smad', 'dlw', 1);
    sigSmcnDLWs = (smcnDLWs - nanmean(smcnDLWs(:))) / nanstd(smcnDLWs(:),1);
    sigSmadDLWs = (smadDLWs - nanmean(smadDLWs(:))) / nanstd(smadDLWs(:),1);

    % check relation between Zi vs signal mean diff, and Zij vs signal amplitude
    checkRelationSubDLWandSignals(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, 'cn');

%    [smrccnDLs, meanSmrccnDL, ~] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'dlcm', 1);
%    [smrccnDLWs, meanSmrccnDLW, ~] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'dlw', 1);
%    [smrcadDLs, meanSmrcadDL, ~] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'dlcm', 1);
%    [smrcadDLWs, meanSmrcadDLW, ~] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'dlw', 1);

    [smcn2DLs, meanSmcn2DL, ~] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'dlcm', 1);
    [smcn2DLWs, meanSmcn2DLW, ~, smcn2SubDLWs] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'dlw', 1);
    [smad2DLs, meanSmad2DL, ~] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'dlcm', 1);
    [smad2DLWs, meanSmad2DLW, ~, smcn2SubDLWs] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'dlw', 1);

%    [smcn3DLs, meanSmcn3DL, ~] = calculateConnectivity(smcn3Signals, roiNames, 'smcn3', 'dlcm', 1);
%    [smcn3DLWs, meanSmcn3DLW, ~] = calculateConnectivity(smcn3Signals, roiNames, 'smcn3', 'dlw', 1);
%    [smad3DLs, meanSmad3DL, ~] = calculateConnectivity(smad3Signals, roiNames, 'smad3', 'dlcm', 1);
%    [smad3DLWs, meanSmad3DLW, ~] = calculateConnectivity(smad3Signals, roiNames, 'smad3', 'dlw', 1);

%    [smcn4DLs, meanSmcn4DL, ~] = calculateConnectivity(smcn4Signals, roiNames, 'smcn4', 'dlcm', 1);
%    [smcn4DLWs, meanSmcn4DLW, ~] = calculateConnectivity(smcn4Signals, roiNames, 'smcn4', 'dlw', 1);
%    [smad4DLs, meanSmad4DL, ~] = calculateConnectivity(smad4Signals, roiNames, 'smad4', 'dlcm', 1);
%    [smad4DLWs, meanSmad4DLW, ~] = calculateConnectivity(smad4Signals, roiNames, 'smad4', 'dlw', 1);

%{
    figure; cnsmcnDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcnDLW);
    figure; adsmadDLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmadDLW);
    figure; cnsmrccnDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmrccnDLW);
    figure; adsmrcadDLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmrcadDLW);
    figure; cnsmcn2DLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn2DLW);
    figure; adsmad2DLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmad2DLW);
    
    m1 = nanmean(meanCnDLW(:));
    m2 = nanmean(meanSmcnDLW(:));
    s1 = nanstd(meanCnDLW(:));
    s2 = nanstd(meanSmcnDLW(:));
    sigCnDLW = (meanCnDLW - m1) ./ s1;
    sigSmcnDLW = (meanSmcnDLW - m2) ./ s2;
    figure; cnsmcnDLWr = plotTwoSignalsCorrelation(sigCnDLW, sigSmcnDLW);
%}    
%    figure; cnsmcn3DLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn3DLW);
%    figure; adsmad3DLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmad3DLW);
%    figure; cnsmcn4DLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn4DLW);
%    figure; adsmad4DLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmad4DLW);

    % check FC of simulated CN and AD
    [cnFCs, meanCnFC, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc', 1);
    [adFCs, meanAdFC, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc', 1);
    [smcnFCs, meanSmcnFC, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'fc', 1);
    [smadFCs, meanSmadFC, ~] = calculateConnectivity(smadSignals, roiNames, 'smad', 'fc', 1);
%    [smrccnFCs, meanSmrccnFC, ~] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'fc', 1);
%    [smrcadFCs, meanSmrcadFC, ~] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'fc', 1);
    [smcn2FCs, meanSmcn2FC, ~] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'fc', 1);
    [smad2FCs, meanSmad2FC, ~] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'fc', 1);
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
    algNum = 6;
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    cosSim(2) = getCosSimilarity(meanCnDLW, meanSmcn2DLW);
    cosSim(3) = getCosSimilarity(meanCnDL, meanSmcnDL);
    cosSim(4) = getCosSimilarity(meanCnDL, meanSmcn2DL);
    cosSim(5) = getCosSimilarity(meanCnFC, meanSmcnFC);
    cosSim(6) = getCosSimilarity(meanCnFC, meanSmcn2FC);
    X = categorical({'dlec-cn-smcn','dlec-cn-smcn2','dlgc-cn-smcn','dlgc-cn-smcn2','fc-cn-smcn','fc-cn-smcn2'});
    figure; bar(X, cosSim); title('cos similarity between mean CN matrix and SimCN by each algorithm');

    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanAdDLW, meanSmadDLW);
    cosSim(2) = getCosSimilarity(meanAdDLW, meanSmad2DLW);
    cosSim(3) = getCosSimilarity(meanAdDL, meanSmadDL);
    cosSim(4) = getCosSimilarity(meanAdDL, meanSmad2DL);
    cosSim(5) = getCosSimilarity(meanAdFC, meanSmadFC);
    cosSim(6) = getCosSimilarity(meanAdFC, meanSmad2FC);
    X = categorical({'dlec-ad-smad','dlec-ad-smad2','dlgc-ad-smad','dlgc-ad-smad2','fc-ad-smad','fc-ad-smad2'});
    figure; bar(X, cosSim); title('cos similarity between mean AD matrix and SimAD by each algorithm');
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
    cosSims = nan(cnSbjNum,9);
    for i=1:cnSbjNum
        cosSims(i,1) = getCosSimilarity(meanCnDLW, cnDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(meanCnDLW, smcnDLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(meanCnDLW, smcn2DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(meanCnDL, cnDLs(:,:,i));
        cosSims(i,5) = getCosSimilarity(meanCnDL, smcnDLs(:,:,i));
        cosSims(i,6) = getCosSimilarity(meanCnDL, smcn2DLs(:,:,i));
        cosSims(i,7) = getCosSimilarity(meanCnFC, cnFCs(:,:,i));
        cosSims(i,8) = getCosSimilarity(meanCnFC, smcnFCs(:,:,i));
        cosSims(i,9) = getCosSimilarity(meanCnFC, smcn2FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:9));
    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
    [sortCosCnDLW,idxCosCnDLW] = sort(cosSims(:,1),'descend');
    [sortCosCnFC,idxCosCnFC] = sort(cosSims(:,7),'descend');

    cosSims = nan(adSbjNum,9);
    for i=1:adSbjNum
        cosSims(i,1) = getCosSimilarity(meanAdDLW, adDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(meanAdDLW, smadDLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(meanAdDLW, smad2DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(meanAdDL, adDLs(:,:,i));
        cosSims(i,5) = getCosSimilarity(meanAdDL, smadDLs(:,:,i));
        cosSims(i,6) = getCosSimilarity(meanAdDL, smad2DLs(:,:,i));
        cosSims(i,7) = getCosSimilarity(meanAdFC, adFCs(:,:,i));
        cosSims(i,8) = getCosSimilarity(meanAdFC, smadFCs(:,:,i));
        cosSims(i,9) = getCosSimilarity(meanAdFC, smad2FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:9));
    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
    [sortCosAdDLW,idxCosAdDLW] = sort(cosSims(:,1),'descend');
    [sortCosAdFC,idxCosAdFC] = sort(cosSims(:,7),'descend');

    % plot correlation between original EC matrix and simulated signals EC matrix
    R = 4; %nodeNum;
    r1 = zeros(R);
    r2 = zeros(R,cnSbjNum);
    for i=1:R
        r1(i) = corr2(squeeze(cnSubDLWs(i,1,:)), squeeze(smcnSubDLWs(i,1,:)));
        figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: cn row=' num2str(i)]);
        for j=1:cnSbjNum
            plotTwoSignalsCorrelation(cnSubDLWs(i,1,j), smcnSubDLWs(i,1,j), [0.1*mod(j,10) 0.2*ceil(j/10) 0.5], 'd', 8);
            plotTwoSignalsCorrelation(cnSubDLWs(i,2:nodeNum+1,j), smcnSubDLWs(i,2:nodeNum+1,j), [0.1*mod(j,10) 0.2*ceil(j/10) 0.8]);
            X = cnSubDLWs(i,2:nodeNum+1,j);
            Y = smcnSubDLWs(i,2:nodeNum+1,j);
            r2(i,j) = corr2(X(~isnan(X)), Y(~isnan(Y)));
        end; hold off;
    end

    % plot box-and-whisker plot of cos similarity between original ec matrix and simulated signals ec matrix
    cosSims = nan(cnSbjNum,6);
    for i=1:cnSbjNum
        cosSims(i,1) = getCosSimilarity(cnDLWs(:,:,i), smcnDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(cnDLWs(:,:,i), smcn2DLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(cnDLs(:,:,i), smcnDLs(:,:,i));
        cosSims(i,4) = getCosSimilarity(cnDLs(:,:,i), smcn2DLs(:,:,i));
        cosSims(i,5) = getCosSimilarity(cnFCs(:,:,i), smcnFCs(:,:,i));
        cosSims(i,6) = getCosSimilarity(cnFCs(:,:,i), smcn2FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    cosSims = nan(cnSbjNum,6);
    for i=1:adSbjNum
        cosSims(i,1) = getCosSimilarity(adDLWs(:,:,i), smadDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(adDLWs(:,:,i), smad2DLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(adDLs(:,:,i), smadDLs(:,:,i));
        cosSims(i,4) = getCosSimilarity(adDLs(:,:,i), smadDLs(:,:,i));
        cosSims(i,5) = getCosSimilarity(adFCs(:,:,i), smadFCs(:,:,i));
        cosSims(i,6) = getCosSimilarity(adFCs(:,:,i), smad2FCs(:,:,i));
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
    [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs, smcnDLs, roiNames, 'cn', 'smcn', 'dlcm');
    [adsmadDLsUt, adsmadDLsUtP, adsmadDLsUtP2] = calculateAlzWilcoxonTest(adDLs, smadDLs, roiNames, 'ad', 'smad', 'dlcm');
    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smcnDLWs, roiNames, 'cn', 'smcn', 'dlw');
    [adsmadDLWsUt, adsmadDLWsUtP, adsmadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smadDLWs, roiNames, 'ad', 'smad', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigCnDLWs, sigSmcnDLWs, roiNames, 'sigcn', 'sigsmcn', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigAdDLWs, sigSmadDLWs, roiNames, 'sigad', 'sigsmad', 'dlw');
%    [cnsmrccnDLWsUt, cnsmrccnDLWsUtP, cnsmrccnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smrccnDLWs, roiNames, 'cn', 'smrccn', 'dlw');
%    [adsmrcadDLWsUt, adsmrcadDLWsUtP, adsmrcadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smrcadDLWs, roiNames, 'ad', 'smrcad', 'dlw');
end

% ==================================================================================================================
function checkRelationSubDLWandSignals(rawSignals, DLWs, subDLWs, simSignals, simDLWs, simSubDLWs, group)
    nodeNum = size(rawSignals{1},1);
    sigLen = size(rawSignals{1},2);
    sbjNum = length(rawSignals);
    R = 8;
    nMax = 20;
    sbjMax = 4;

    % checking signal parallel shift effect for Zi, Zij and ECij'
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);
        smSi = simSignals{k};
        smEC = simDLWs(:,:,k);
        smSubEC = simSubDLWs(:,:,k);
        
        outfName = ['results/adsim2-checkRelation-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 20;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end
    
            Zi2 = zeros(R, nMax, 17);
            X = zeros(R, nMax, 17);
            Zij2 = zeros(R, nMax, 17, nodeNum);
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            [siOrg, sig, c, maxsi, minsi] = convert2SigmoidSignal(rawSignals{k});

            % training options for DLCM network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 3);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.05, ...
                'Verbose',false);
                
            for i=1:R 
                Zi(i) = subEC(i,1); % original Zi value
                Zij(i,:) = subEC(i,2:end); % original Zij value
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.inControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                for a=0:16
                    dx = (-0.32 + 0.04*a);
                    si = siOrg;
                    si(i,:) = si(i,:) + dx;

                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); f.inSignal(:,1:end-1)];
                    filter = repmat(f.inControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    subEC2 = cell(nMax,1);
%                    for n=1:nMax % traial
                    parfor n=1:nMax % traial
                        netDLCM = initDlcmNetwork(si, f.inSignal, [], f.inControl); 

                        disp(['training ' num2str(k) '-' num2str(i) ' dx=' num2str(dx) ' n:' num2str(n)]);
                        [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                        % predict DLCM network
                        subEC2{n} = predict(nodeNetwork, Si1);
                    end
                    for n=1:nMax
                        Zi2(i, n, a+1) = subEC2{n}(1);
                        Zij2(i, n, a+1,:) = subEC2{n}(2:end);
                        X(i, n, a+1) = dx;
                    end
                end
            end
            save(outfName, 'Zi', 'Zi2', 'Zij2', 'X');

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
        
        % plot result -- Zi2 vs dx
        figure;
        for i=1:R
            x=X(i,:,:);
            y=Zi2(i,:,:);
            hold on; scatter(x(:),y(:),3); hold off;
        end
        daspect([1 1 1]); title(['sbj' num2str(k) ' Zi vs dx']);

        % plot result -- Zij2(1:16) vs dx
        figure;
        for i=1:1
            for j=1:16
                x=X(i,:,:);
                y=Zij2(i,:,:,j);
                hold on; scatter(x(:),y(:),3); hold off;
            end
        end
        daspect([1 1 1]); title(['sbj' num2str(k) ' Zij vs dx']);

        % plot result -- Zi - Zij2(1:16) vs dx
        figure;
        for i=1:1
            for j=1:16
                x=X(i,:,:);
                y=Zi2(i,:,:) - Zij2(i,:,:,j);
                hold on; scatter(x(:),y(:),3); hold off;
            end
        end
        daspect([1 1 1]); title(['sbj' num2str(k) ' (Zi - Zij) vs dx']);
    end
    
    % checking signal amplitude change effect for Zi, Zij and ECij'
    amps = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5, 2, 3, 5, 8];
    ampsLen = length(amps);
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);

        outfName = ['results/adsim2-checkRelation2-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 20;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end
    
            Zi2 = zeros(R, nMax, ampsLen);
            X = zeros(R, nMax, ampsLen);
            Zij2 = zeros(R, nMax, ampsLen, nodeNum);
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            [siOrg, sig, c, maxsi, minsi] = convert2SigmoidSignal(rawSignals{k});

            % training options for DLCM network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 3);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.05, ...
                'Verbose',false);
                
            for i=1:R 
                Zi(i) = subEC(i,1); % original Zi value
                Zij(i,:) = subEC(i,2:end); % original Zij value
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.inControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                for a=1:ampsLen
                    si = siOrg;
                    amp = amps(a);
                    m = nanmean(si(i,:));
                    si(i,:) = (si(i,:)-m) .* amp + m;

                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); f.inSignal(:,1:end-1)];
                    filter = repmat(f.inControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    subEC2 = cell(nMax,1);
%                    for n=1:nMax % traial
                    parfor n=1:nMax % traial
                        netDLCM = initDlcmNetwork(si, f.inSignal, [], f.inControl); 

                        disp(['training ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);
                        [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                        % predict DLCM network
                        subEC2{n} = predict(nodeNetwork, Si1);
                    end
                    for n=1:nMax
                        Zi2(i, n, a) = subEC2{n}(1);
                        Zij2(i, n, a,:) = subEC2{n}(2:end);
                        X(i, n, a) = amp;
                    end
                end
            end
            save(outfName, 'Zi', 'Zi2', 'Zij2', 'X');

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
        
        % plot result -- Zi2 vs dx
        figure;
        for i=1:R
            x=X(i,:,:);
            y=Zi2(i,:,:);
            hold on; scatter(x(:),y(:),3); hold off;
        end
        daspect([1 1 1]); title(['sbj' num2str(k) ' Zi vs dx']);

        % plot result -- Zij2(1:16) vs dx
        figure;
        for i=1:1
            for j=1:16
                x=X(i,:,:);
                y=Zij2(i,:,:,j);
                hold on; scatter(x(:),y(:),3); hold off;
            end
        end
        daspect([1 1 1]); title(['sbj' num2str(k) ' Zij vs dx']);

        % plot result -- Zi - Zij2(1:16) vs dx
        figure;
        for i=1:1
            for j=1:16
                x=X(i,:,:);
                y=Zi2(i,:,:) - Zij2(i,:,:,j);
                hold on; scatter(x(:),y(:),3); hold off;
            end
        end
        daspect([1 1 1]); title(['sbj' num2str(k) ' (Zi - Zij) vs dx']);
    end
end

function [smDLWs, bSignals] = retrainDLCMAndECmultiPattern(rawSignals, DLWs, subDLWs, S2, IS2, roiNames, group)
    nodeNum = size(rawSignals{1},1);
    sigLen = size(rawSignals{1},2);
    sbjNum = length(rawSignals);
    nanx = eye(nodeNum);
    nanx(nanx==1) = NaN;

    R = nodeNum;
    JMAX = 7;
    k1 = floor(101/20)+1;
    r1 = zeros(JMAX+1,k1,R);
    r2 = zeros(JMAX+1,k1,R,sbjNum);
    r3 = zeros(JMAX+1,k1,sbjNum);
    h1 = zeros(JMAX+1,k1,nodeNum,nodeNum);
    p1 = zeros(JMAX+1,k1,nodeNum,nodeNum);
    h1c = zeros(JMAX+1,k1);
    r1b = zeros(JMAX+1,k1,R);
    r2b = zeros(JMAX+1,k1,R,sbjNum);
    r3b = zeros(JMAX+1,k1,sbjNum);
    h1b = zeros(JMAX+1,k1,nodeNum,nodeNum);
    p1b = zeros(JMAX+1,k1,nodeNum,nodeNum);
    h1bc = zeros(JMAX+1,k1);

    Zi = repmat(subDLWs(:,1,:),[1 nodeNum 1]);
    Zij = subDLWs(:,2:nodeNum+1,:);
    DLWsR = Zi - Zij;
%    DLWs = abs(DLWsR);

    for ii=0:0
        for exRate=2:2 %0:6
            for k=1:1 %1:20:101
                ECij = Zij - DLWsR * exRate * 0.2;
                teaches = [];
                S3 = [];
                IS3 = [];
                signals = {};
                for i=1:sbjNum
                    [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(rawSignals{i});
                    signals{end+1} = si;
                    teach = [repmat(subDLWs(:,1,i),[1 k 1]) ECij(:,:,i)]; % indivisual part of teaching data
                    teaches(:,:,i) = [teach si(:,2:end)];
                    S3(:,:,i) = [repmat(S2(:,1),[1 k]) S2(:,2:nodeNum+1) si(:,1:end-1)];
                    IS3(:,:,i) = [repmat(IS2(:,1),[1 k]) IS2(:,2:nodeNum+1) rand(nodeNum,sigLen-1)];
                end
                
                name = [group '-' num2str(ii) '-' num2str(exRate) '-' num2str(k) 'ns'];
                bname = [group 'b-' num2str(ii) '-' num2str(exRate) '-' num2str(k) 'ns'];
                ecname = [group '-' num2str(ii) '-' num2str(exRate)  '-' num2str(k) 'ec'];

                % train DLCM network by signals & expanded EC
                [bDLWs, meanbDLWns, stdbDLWns, bSubDLWs] = retrainDLCMAndEC(teaches, S3, IS3, roiNames, name);
                
                k1 = floor(k/20)+1;
                for b=1:R
                    r1b(exRate+1,k1,b) = corr2(squeeze(subDLWs(b,1,:)), squeeze(bSubDLWs(b,1,:)));
%                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' name ' row=' num2str(b)]);
                    for a=1:sbjNum
%                        plotTwoSignalsCorrelation(subDLWs(b,1,a), bSubDLWs(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
%                        plotTwoSignalsCorrelation(subDLWs(b,2:77,a), bSubDLWs(b,2:77,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]);
                        r2b(exRate+1,k1,b,a) = corr2(subDLWs(b,2:1+nodeNum,a), bSubDLWs(b,2:1+nodeNum,a));
                    end; hold off;
                end
%                figure; hold on; plot([0 0.5], [0 0.5],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' ecname ' row=' num2str(b)]);
                for a=1:sbjNum
                    X = DLWs(1:R,1:R,a)+nanx(1:R,1:R);
                    Y = bDLWs(1:R,1:R,a);
%                    plotTwoSignalsCorrelation(X, Y, [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                    r3b(exRate+1,k1,a) = corr2(X(~isnan(X(:))), Y(~isnan(Y(:))));
                end; hold off;
%                calculateAlzWilcoxonTest(subDLWs, bSubDLWs, roiNames, 'ns', name, 'dlw', 1, 'ranksum');
                [h1b(exRate+1,k1,:,:), p1b(exRate+1,k1,:,:), ~] = calculateAlzWilcoxonTest(DLWs, bDLWs, roiNames, 'ec', ecname, 'dlw', 1, 'ranksum', 0);
                h1bc(exRate+1,k1) = length(find(h1b(exRate+1,k1,1:R,:)>0));

                % simulate from first frame with the DLCM network trained by signals & expanded EC
                [bDLWs, bSignals, bSubDLWs] = simulateNodeSignals(signals, roiNames, name, 'dlw', name, 1, [k+nodeNum+1:k+nodeNum+sigLen-1]);
                [smDLs, meanSmDL, ~] = calculateConnectivity(bSignals, roiNames, bname, 'dlcm', 1);
                [smDLWs, meanSmDLW, ~, smSubDLWs] = calculateConnectivity(bSignals, roiNames, bname, 'dlw', 1);

                for b=1:R
                    r1(exRate+1,k1,b) = corr2(squeeze(subDLWs(b,1,:)), squeeze(smSubDLWs(b,1,:)));
%                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' name ' row=' num2str(b)]);
                    for a=1:sbjNum
%                        plotTwoSignalsCorrelation(subDLWs(b,1,a), smSubDLWs(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
%                        plotTwoSignalsCorrelation(subDLWs(b,2:77,a), smSubDLWs(b,2:77,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]);
                        r2(exRate+1,k1,b,a) = corr2(subDLWs(b,2:1+nodeNum,a), smSubDLWs(b,2:1+nodeNum,a));
                    end; hold off;
                end
%                figure; hold on; plot([0 0.5], [0 0.5],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' ecname ' row=' num2str(b)]);
                for a=1:sbjNum
                    X = DLWs(1:R,1:R,a)+nanx(1:R,1:R);
                    Y = smDLWs(1:R,1:R,a);
%                    plotTwoSignalsCorrelation(X, Y, [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                    r3(exRate+1,k1,a) = corr2(X(~isnan(X(:))), Y(~isnan(Y(:))));
                end; hold off;
%                calculateAlzWilcoxonTest(subDLWs, smSubDLWs, roiNames, 'ns', name, 'dlw', 1, 'ranksum');
                [h1(exRate+1,k1,:,:), p1(exRate+1,k1,:,:), ~] = calculateAlzWilcoxonTest(DLWs, smDLWs, roiNames, 'ec', ecname, 'dlw', 1, 'ranksum', 0);
                h1c(exRate+1,k1) = length(find(h1(exRate+1,k1,1:R,:)>0));
            end
        end
    end
    r1m = nanmean(r1,3);
    r2m = nanmean(nanmean(r2,4),3);
    r3m = nanmean(r3,3);
    p1m = nanmean(nanmean(p1(:,:,1:R,:),4),3);
    r1bm = nanmean(r1b,3);
    r2bm = nanmean(nanmean(r2b,4),3);
    r3bm = nanmean(r3b,3);
    p1bm = nanmean(nanmean(p1b(:,:,1:R,:),4),3);
end

function [weights, meanWeights, stdWeights, subweights] = retrainDLCMAndEC(teachSignals, nodeSignals, exSignals, roiNames, group)
    ROWNUM = size(teachSignals,1);
    COLNUM = size(teachSignals,2);
    sbjNum = size(teachSignals,3);
    weights = zeros(ROWNUM, ROWNUM, sbjNum);

    outfName = ['results/adsim2-retrain-' group '-roi' num2str(ROWNUM) '.mat'];
    if exist(outfName, 'file')
        f=load(outfName);
        weights = f.weights;
        subweights = f.subweights;
        meanWeights = nanmean(f.weights, 3);
        stdWeights = nanstd(f.weights, 1, 3);
        return;
    end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 14;

    if NumProcessors > 1
        try
            disp('Destroing any existance matlab pool session');
            parpool('close');
        catch
            disp('No matlab pool session found');
        end
        parpool(NumProcessors);
    end

    % init params
    sigLen = size(nodeSignals,2);
    inControl = eye(ROWNUM);

%    for i=1:sbjNum
    parfor i=1:sbjNum
        dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROWNUM) '-net' num2str(i) '.mat'];
        if exist(dlcmName, 'file')
            f=load(dlcmName);
            netDLCM = f.netDLCM;
        else
            if size(nodeSignals,3) > 1
                si = nodeSignals(:,:,i);
                inSignal = exSignals(:,:,i);
            else
                si = nodeSignals;
                inSignal = exSignals;
            end
            % init DLCM network
            netDLCM = initDlcmNetwork(si, inSignal, [], inControl);

            % training DLCM network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 2);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.05, ...
                'Verbose',false);
        %            'Plots','training-progress');

            disp('start training');
            for j=1:ROWNUM
                nodeTeach = teachSignals(j,1:end,i);
                nodeInput = [si; inSignal];
                if ~isempty(inControl)
                    filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(ROWNUM+1:end,:) = nodeInput(ROWNUM+1:end,:) .* filter;
                end
                idx = find(isnan(nodeTeach));
                nodeTeach(:,idx) = [];
                nodeInput(:,idx) = [];
                [netDLCM.nodeNetwork{j}, netDLCM.trainInfo{j}] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{j}, options);
                disp(['virtual alzheimer (' group ') training node ' num2str(i) '-' num2str(j) ' rmse=' num2str(netDLCM.trainInfo{j}.TrainingRMSE(maxEpochs))]);
            end

            parsavedlsm(dlcmName, netDLCM, si, inSignal, inControl, options);
        end

        % recalculate EC
        [weights(:,:,i), subweights(:,:,i)] = calcDlcmEC(netDLCM, [], inControl);
    end
    save(outfName, 'weights', 'roiNames', 'subweights');
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
end

function parsavedlsm(dlcmName, netDLCM, si, inSignal, inControl, options)
    save(dlcmName, 'netDLCM', 'si', 'inSignal', 'inControl', 'options');
end

function out = expandAmplitude(signals, rate)
    m = nanmean(signals(:));
    d = signals - m;
    out = d * rate + m;
end

function out = expandAmplitude2(signals, rate)
    m = nanmean(signals,2);
    d = signals - m;
    out = d * rate + m;
end

function [ECs, simSignals, subECs] = simulateNodeSignals(signals, roiNames, group, algorithm, orgGroup, isRaw, inSiRange)
    if nargin < 7
        inSiRange = 0;
    end
    if nargin < 6
        isRaw = 0;
    end
    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    outfName = ['results/adsim2-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 14;

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
        case {'dlw','dlwrc'}
            if strcmp(algorithm, 'dlw')
                dlcmName = ['results/ad-dlcm-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            else
                dlcmName = ['results/ad-dlcmrc-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            end
            f = load(dlcmName);
            if isRaw
                si = signals{i};
                sig=0; c=0; maxsi=0; minsi=0;
            else
                [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
            end
            if inSiRange > 0
                inSignal = f.inSignal(:,inSiRange);
            else
                inSignal = f.inSignal;
            end
            [Y, time] = simulateDlcmNetwork(si, inSignal, [], f.inControl, f.netDLCM);
            [ec, subECs(:,:,i)] = calcDlcmEC(f.netDLCM, [], f.inControl);
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
