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

    % get DLCM-EC, DLCM-GC & FC from original CN and AD (normal and recovery training)
    for j=1:1
        [cnDLs{j}, meanCnDL{j}, stdCnDL{j}] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 0, j, 1);
        [adDLs{j}, meanAdDL{j}, stdAdDL{j}] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 0, j, 1);
        sigCnDLs{j} = (cnDLs{j} - nanmean(cnDLs{j}(:))) / nanstd(cnDLs{j}(:),1);
        sigAdDLs{j} = (adDLs{j} - nanmean(adDLs{j}(:))) / nanstd(adDLs{j}(:),1);

        [cnFCs{j}, meanCnFC{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc', 1);
        [adFCs{j}, meanAdFC{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc', 1);
    end

    % --------------------------------------------------------------------------------------------------------------
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
        [cnDLWs{j}, smcnSignals{j}, cnSubDLWs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn', 0, 0, i, 1);
        [adDLWs{j}, smadSignals{j}, adSubDLWs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad', 0, 0, i, 1);
        sigCnDLWs{j} = (cnDLWs{j} - nanmean(cnDLWs{j}(:))) / nanstd(cnDLWs{j}(:),1);
        sigAdDLWs{j} = (adDLWs{j} - nanmean(adDLWs{j}(:))) / nanstd(adDLWs{j}(:),1);
        meanCnDLW{j} = nanmean(cnDLWs{j},3);
        meanAdDLW{j} = nanmean(adDLWs{j},3);

        % DLCM(j) linear auto exogenous 
        [cnDLW2s{j}, sm2cnSignals{j}, cnSubDLW2s{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn', 0, 0, i, 1, []);
        [adDLW2s{j}, sm2adSignals{j}, adSubDLW2s{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad', 0, 0, i, 1, []);
        sigCnDLW2s{j} = (cnDLW2s{j} - nanmean(cnDLW2s{j}(:))) / nanstd(cnDLW2s{j}(:),1);
        sigAdDLW2s{j} = (adDLW2s{j} - nanmean(adDLW2s{j}(:))) / nanstd(adDLW2s{j}(:),1);
        meanCnDLW2{j} = nanmean(cnDLW2s{j},3);
        meanAdDLW2{j} = nanmean(adDLW2s{j},3);

        % mvar(j) auto exogenous 
        [cnMVARECs{j}, smmvcnSignals{j}, cnSubMVARECs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'mvarec', 'cn', 0, 0, i, 1);
        [adMVARECs{j}, smmvadSignals{j}, adSubMVARECs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'mvarec', 'ad', 0, 0, i, 1);
        sigCnDLWs{j} = (cnMVARECs{j} - nanmean(cnMVARECs{j}(:))) / nanstd(cnMVARECs{j}(:),1);
        sigAdDLWs{j} = (adMVARECs{j} - nanmean(adMVARECs{j}(:))) / nanstd(adMVARECs{j}(:),1);
        meanCnDLW{j} = nanmean(cnMVARECs{j},3);
        meanAdDLW{j} = nanmean(adMVARECs{j},3);
    end

    % --------------------------------------------------------------------------------------------------------------
    % check DLCM-EC and DLCM-GC of simulated CN and AD
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [smcnDLs{j}, meanSmcnDL{j}, ~] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'dlcm', 1, j, 0);
        [smcnDLWs{j}, meanSmcnDLW{j}, ~, smcnSubDLWs{j}] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'dlw', 1, j, 0);
        % DLCM(j) linear no exogenous 
        [sm2cnDLs{j}, meanSm2cnDL{j}, ~] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'dlcm', 1, j, 0, []);
        [sm2cnDLWs{j}, meanSm2cnDLW{j}, ~, sm2cnSubDLWs{j}] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'dlw', 1, j, 0, []);
        % mvar(j) no exogenous 
        [smcnDL2s{j}, meanSmcnDL2{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smcn', 'dlcm', 1, j, 0);
        [smcnDLW2s{j}, meanSmcnDLW2{j}, ~, smcnSubDLW2s{j}] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smcn', 'dlw', 1, j, 0);
    end
    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [smcnDLs{j}, meanSmcnDL{j}, ~] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'dlcm', 1, i, 1);
        [smcnDLWs{j}, meanSmcnDLW{j}, ~, smcnSubDLWs{j}] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'dlw', 1, i, 1);
        % DLCM(j) linear auto exogenous 
        [sm2cnDLs{j}, meanSm2cnDL{j}, ~] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'dlcm', 1, i, 1, []);
        [sm2cnDLWs{j}, meanSm2cnDLW{j}, ~, sm2cnSubDLWs{j}] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'dlw', 1, i, 1, []);
        % mvar(j) auto exogenous 
        [smcnDL2s{j}, meanSmcnDL2{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smcn', 'dlcm', 1, i, 1);
        [smcnDLW2s{j}, meanSmcnDLW2{j}, ~, smcnSubDLW2s{j}] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smcn', 'dlw', 1, i, 1);
    end

    % check FC of simulated CN and AD
    for j=1:maxLag
        % DLCM(j) auto exogenous 
        [smcnFCs{j}, meanSmcnFC{j}, ~] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'fc', 1, j, 0);
        [smadFCs{j}, meanSmadFC{j}, ~] = calculateConnectivity(smadSignals{j}, roiNames, 'smad', 'fc', 1, j, 0);
        % DLCM(j) auto exogenous 
        [sm2cnFCs{j}, meanSm2cnFC{j}, ~] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'fc', 1, j, 0, []);
        [sm2adFCs{j}, meanSm2adFC{j}, ~] = calculateConnectivity(sm2adSignals{j}, roiNames, 'smad', 'fc', 1, j, 0, []);
        % mvar(j) auto exogenous 
        [smmvcnFCs{j}, meanSmmvcnFC{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smcn', 'fc', 1, j, 0);
        [smmvadFCs{j}, meanSmmvadFC{j}, ~] = calculateConnectivity(smmvadSignals{j}, roiNames, 'smad', 'fc', 1, j, 0);
    end
    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [smcnFCs{j}, meanSmcnFC{j}, ~] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'fc', 1, i, 1);
        [smadFCs{j}, meanSmadFC{j}, ~] = calculateConnectivity(smadSignals{j}, roiNames, 'smad', 'fc', 1, i, 1);
        % DLCM(j) auto exogenous 
        [sm2cnFCs{j}, meanSm2cnFC{j}, ~] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'fc', 1, i, 1, []);
        [sm2adFCs{j}, meanSm2adFC{j}, ~] = calculateConnectivity(sm2adSignals{j}, roiNames, 'smad', 'fc', 1, i, 1, []);
        % mvar(j) auto exogenous 
        [smmvcnFCs{j}, meanSmmvcnFC{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smcn', 'fc', 1, i, 1);
        [smmvadFCs{j}, meanSmmvadFC{j}, ~] = calculateConnectivity(smmvadSignals{j}, roiNames, 'smad', 'fc', 1, i, 1);
    end

    % --------------------------------------------------------------------------------------------------------------
    % plot correlation and cos similarity
    cosSim = zeros(60,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    cosSim(2) = getCosSimilarity(meanCnDLW, meanSmcn2DLW);
    cosSim(11) = getCosSimilarity(meanCnDL, meanSmcnDL);
    cosSim(12) = getCosSimilarity(meanCnDL, meanSmcn2DL);
    cosSim(21) = getCosSimilarity(meanCnFC, meanSmcnFC);
    cosSim(22) = getCosSimilarity(meanCnFC, meanSmcn2FC);
    figure; bar(cosSim); title('cos similarity between mean CN matrix and SimCN by each algorithm');

    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanAdDLW, meanSmadDLW);
    cosSim(2) = getCosSimilarity(meanAdDLW, meanSmad2DLW);
    cosSim(11) = getCosSimilarity(meanAdDL, meanSmadDL);
    cosSim(12) = getCosSimilarity(meanAdDL, meanSmad2DL);
    cosSim(21) = getCosSimilarity(meanAdFC, meanSmadFC);
    cosSim(22) = getCosSimilarity(meanAdFC, meanSmad2FC);
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
        cosSims(i,11) = getCosSimilarity(meanCnDL, cnDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(meanCnDL, smcnDLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(meanCnFC, cnFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(meanCnFC, smcnFCs(:,:,i));
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
        cosSims(i,11) = getCosSimilarity(meanAdDL, adDLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(meanAdFC, adFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(meanAdFC, smadFCs(:,:,i));
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
        cosSims(i,11) = getCosSimilarity(cnDLs(:,:,i), smcnDLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(cnFCs(:,:,i), smcnFCs(:,:,i));
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
        cosSims(i,11) = getCosSimilarity(adDLs(:,:,i), smadDLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(adFCs(:,:,i), smadFCs(:,:,i));
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
    [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs, smcnDLs, roiNames, 'cn', 'smcn', 'dlcm');
    [adsmadDLsUt, adsmadDLsUtP, adsmadDLsUtP2] = calculateAlzWilcoxonTest(adDLs, smadDLs, roiNames, 'ad', 'smad', 'dlcm');
    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smcnDLWs, roiNames, 'cn', 'smcn', 'dlw');
    [adsmadDLWsUt, adsmadDLWsUtP, adsmadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smadDLWs, roiNames, 'ad', 'smad', 'dlw');
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
