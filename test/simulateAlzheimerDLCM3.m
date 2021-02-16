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

        % mpcvar(j) no exogenous 
        [cnMPCVARECs{j}, smmpvcnSignals{j}, cnSubMPCVARECs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'mpcvarec', 'cn', 0, 0, j, 0);
        [adMPCVARECs{j}, smmpvadSignals{j}, adSubMPCVARECs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'mpcvarec', 'ad', 0, 0, j, 0);
        sigCnMPCVARECs{j} = (cnMPCVARECs{j} - nanmean(cnMPCVARECs{j}(:))) / nanstd(cnMPCVARECs{j}(:),1);
        sigAdMPCVARECs{j} = (adMPCVARECs{j} - nanmean(adMPCVARECs{j}(:))) / nanstd(adMPCVARECs{j}(:),1);
        meanCnMPCVAREC{j} = nanmean(cnMPCVARECs{j},3);
        meanAdMPCVAREC{j} = nanmean(adMPCVARECs{j},3);
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

        % mpcvar(j) auto exogenous 
        [cnMPCVARECs{j}, smmpvcnSignals{j}, cnSubMPCVARECs{j}] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'mpcvarec', 'cn', 0, 0, i, 1);
        [adMPCVARECs{j}, smmpvadSignals{j}, adSubMPCVARECs{j}] = simulateNodeSignals(adSignals, roiNames, 'ad', 'mpcvarec', 'ad', 0, 0, i, 1);
        sigCnMPCVARECs{j} = (cnMPCVARECs{j} - nanmean(cnMPCVARECs{j}(:))) / nanstd(cnMPCVARECs{j}(:),1);
        sigAdMPCVARECs{j} = (adMPCVARECs{j} - nanmean(adMPCVARECs{j}(:))) / nanstd(adMPCVARECs{j}(:),1);
        meanCnMPCVAREC{j} = nanmean(cnMPCVARECs{j},3);
        meanAdMPCVAREC{j} = nanmean(adMPCVARECs{j},3);
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
        [smmvcnDLs{j}, meanSmmvcnDL{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smmvcn', 'dlcm', 1, j, 0);
        [smmvcnDLWs{j}, meanSmmvcnDLW{j}, ~, smmvcnSubDLW2s{j}] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smmvcn', 'dlw', 1, j, 0);
        % mvar(j) no exogenous 
        [smmpvcnDLs{j}, meanSmmpvcnDL{j}, ~] = calculateConnectivity(smmpvcnSignals{j}, roiNames, 'smmpvcn', 'dlcm', 1, j, 0);
        [smmpvcnDLWs{j}, meanSmmpvcnDLW{j}, ~, smmpvcnSubDLW2s{j}] = calculateConnectivity(smmpvcnSignals{j}, roiNames, 'smmpvcn', 'dlw', 1, j, 0);
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
        [smmvcnDLs{j}, meanSmmvcnDL{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smmvcn', 'dlcm', 1, i, 1);
        [smmvcnDLWs{j}, meanSmmvcnDLW{j}, ~, smcnSubDLW2s{j}] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smmvcn', 'dlw', 1, i, 1);
        % mvar(j) auto exogenous 
        [smmpvcnDLs{j}, meanSmmpvcnDL{j}, ~] = calculateConnectivity(smmpvcnSignals{j}, roiNames, 'smmpvcn', 'dlcm', 1, i, 1);
        [smmpvcnDLWs{j}, meanSmmpvcnDLW{j}, ~, smmpvcnSubDLW2s{j}] = calculateConnectivity(smmpvcnSignals{j}, roiNames, 'smmpvcn', 'dlw', 1, i, 1);
    end

    % check FC of simulated CN and AD
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [smcnFCs{j}, meanSmcnFC{j}, ~] = calculateConnectivity(smcnSignals{j}, roiNames, 'smcn', 'fc', 1, j, 0);
        [smadFCs{j}, meanSmadFC{j}, ~] = calculateConnectivity(smadSignals{j}, roiNames, 'smad', 'fc', 1, j, 0);
        % DLCM(j) no exogenous 
        [sm2cnFCs{j}, meanSm2cnFC{j}, ~] = calculateConnectivity(sm2cnSignals{j}, roiNames, 'smcn', 'fc', 1, j, 0, []);
        [sm2adFCs{j}, meanSm2adFC{j}, ~] = calculateConnectivity(sm2adSignals{j}, roiNames, 'smad', 'fc', 1, j, 0, []);
        % mvar(j) no exogenous 
        [smmvcnFCs{j}, meanSmmvcnFC{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smmvcn', 'fc', 1, j, 0);
        [smmvadFCs{j}, meanSmmvadFC{j}, ~] = calculateConnectivity(smmvadSignals{j}, roiNames, 'smmvad', 'fc', 1, j, 0);
        % mpcvar(j) no exogenous 
        [smmpvcnFCs{j}, meanSmmpvcnFC{j}, ~] = calculateConnectivity(smmpvcnSignals{j}, roiNames, 'smmpvcn', 'fc', 1, j, 0);
        [smmpvadFCs{j}, meanSmmpvadFC{j}, ~] = calculateConnectivity(smmpvadSignals{j}, roiNames, 'smmpvad', 'fc', 1, j, 0);
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
        [smmvcnFCs{j}, meanSmmvcnFC{j}, ~] = calculateConnectivity(smmvcnSignals{j}, roiNames, 'smmvcn', 'fc', 1, i, 1);
        [smmvadFCs{j}, meanSmmvadFC{j}, ~] = calculateConnectivity(smmvadSignals{j}, roiNames, 'smmvad', 'fc', 1, i, 1);
        % mpcvar(j) auto exogenous 
        [smmpvcnFCs{j}, meanSmmpvcnFC{j}, ~] = calculateConnectivity(smmpvcnSignals{j}, roiNames, 'smmpvcn', 'fc', 1, i, 1);
        [smmpvadFCs{j}, meanSmmpvadFC{j}, ~] = calculateConnectivity(smmpvadSignals{j}, roiNames, 'smmpvad', 'fc', 1, i, 1);
    end

    % --------------------------------------------------------------------------------------------------------------
    % box plot of signal MAEs 
    maes = nan(cnSbjNum,40);
    for j=1:cnSbjNum
        si = convert2SigmoidSignal(cnSignals{j});
        for i=1:maxLag*2
            maes(j,i) = getTwoSignalsError(si, smcnSignals{i}{j}); % non linear (no ex, ex)
            maes(j,10+i) = getTwoSignalsError(si, sm2cnSignals{i}{j}); % linear (no ex, ex)
            maes(j,20+i) = getTwoSignalsError(si, smmvcnSignals{i}{j}); % linear (no ex, ex)
            maes(j,30+i) = getTwoSignalsError(si, smmpvcnSignals{i}{j}); % linear (no ex, ex)
        end
    end
    figure; boxplot(maes); title('MAEs between mean cnSignals and each algorithm signals');

    % plot correlation and cos similarity
    cosSim = zeros(90,1);
    for k=1:maxLag*2
        i=k;
        cosSim(i) = getCosSimilarity(meanCnDLW{1}, meanSmcnDLW{k}); % non linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanCnDLW{1}, meanSm2cnDLW{k}); i=i+10; % linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanCnDLW{1}, meanSmmvcnDLW{k}); i=i+10; % mvar linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanCnDLW{1}, meanSmmpvcnDLW{k}); i=i+10; % mpcvar linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanCnDL{1}, meanSmcnDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnDL{1}, meanSm2cnDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnDL{1}, meanSmmvcnDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnDL{1}, meanSmmpvcnDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnFC{1}, meanSmcnFC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnFC{1}, meanSm2cnFC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnFC{1}, meanSmmvcnFC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanCnFC{1}, meanSmmpvcnFC{k}); i=i+10;
    end
    figure; bar(cosSim); title('cos similarity between mean CN matrix and SimCN by each algorithm');
%{
    cosSim = zeros(algNum,1);
    figure; bar(cosSim); title('cos similarity between mean AD matrix and SimAD by each algorithm');
%}  
    % plot box-and-whisker plot of cos similarity between mean ec matrix and each subject ec
    cosSims = nan(cnSbjNum,90);
    for j=1:cnSbjNum
        for k=1:maxLag*2
            i=k;
            cosSims(j,i) = getCosSimilarity(meanCnDLW{1}, smcnDLWs{k}(:,:,j)); % non linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanCnDLW{1}, sm2cnDLWs{k}(:,:,j)); i=i+10; % linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanCnDLW{1}, smmvcnDLWs{k}(:,:,j)); i=i+10; % mvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanCnDLW{1}, smmpvcnDLWs{k}(:,:,j)); i=i+10; % mpcvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanCnDL{1}, smcnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnDL{1}, sm2cnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnDL{1}, smmvcnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnDL{1}, smmpvcnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnFC{1}, smcnFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnFC{1}, sm2cnFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnFC{1}, smmvcnFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanCnFC{1}, smmpvcnFCs{k}(:,:,j)); i=i+10;
        end
%        cosSims(j,1) = getCosSimilarity(meanCnDLW, cnDLWs(:,:,i));
%        cosSims(i,11) = getCosSimilarity(meanCnDL, cnDLs(:,:,i));
%        cosSims(i,21) = getCosSimilarity(meanCnFC, cnFCs(:,:,i));
    end
    figure; boxplot(cosSims); title('cos similarity between mean CN matrix and each subject EC by each algorithm');
%    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:30));
%    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
%    [sortCosCnDLW,idxCosCnDLW] = sort(cosSims(:,1),'descend');
%    [sortCosCnFC,idxCosCnFC] = sort(cosSims(:,21),'descend');

    cosSims = nan(adSbjNum,90);
    for j=1:adSbjNum
        for i=1:maxLag*2
        end
    end
%    figure; boxplot(cosSims);

    % plot correlation between original EC matrix and simulated signals EC matrix
%{
    for i=1:1 %cnSbjNum
        % calc & plot correlation of original Zi vs simulating Zi
        [corrZi, corrZij] = calcCorrelationZiZij(cnSubDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum);
        plotCorrelationZiZij(cnDLWs(:,:,i), cnSubDLWs(:,:,i), smcnDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum, ['sbj' num2str(i)], 'original', 'simulating');
    end
%}
    % plot box-and-whisker plot of cos similarity between original ec matrix and simulated signals ec matrix
    cosSims = nan(cnSbjNum,90);
    for j=1:cnSbjNum
        for k=1:maxLag*2
            i=k;
            cosSims(j,i) = getCosSimilarity(cnDLWs{1}(:,:,j), smcnDLWs{k}(:,:,j)); % non linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(cnDLWs{1}(:,:,j), sm2cnDLWs{k}(:,:,j)); i=i+10; % linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(cnDLWs{1}(:,:,j), smmvcnDLWs{k}(:,:,j)); i=i+10; % mvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(cnDLWs{1}(:,:,j), smmpvcnDLWs{k}(:,:,j)); i=i+10; % mpcvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(cnDLs{1}(:,:,j), smcnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnDLs{1}(:,:,j), sm2cnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnDLs{1}(:,:,j), smmvcnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnDLs{1}(:,:,j), smmpvcnDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnFCs{1}(:,:,j), smcnFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnFCs{1}(:,:,j), sm2cnFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnFCs{1}(:,:,j), smmvcnFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(cnFCs{1}(:,:,j), smmpvcnFCs{k}(:,:,j)); i=i+10;
        end
    end
    figure; boxplot(cosSims); title('cos similarity between mean CN matrix and each subject EC by each algorithm');
%{
    cosSims = nan(cnSbjNum,30);
    for j=1:adSbjNum
        for i=1:maxLag*2
            cosSims(i,1) = getCosSimilarity(adDLWs(:,:,i), smadDLWs(:,:,i));
            cosSims(i,11) = getCosSimilarity(adDLs(:,:,i), smadDLs(:,:,i));
            cosSims(i,21) = getCosSimilarity(adFCs(:,:,i), smadFCs(:,:,i));
        end
    end
    figure; boxplot(cosSims);
%}
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
    for i=1:3 %maxLag*2
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs{1}, smcnFCs{i}, roiNames, 'cn', ['sm' num2str(i) 'cn'], 'fc');
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs{1}, sm2cnFCs{i}, roiNames, 'cn', ['smlin' num2str(i) 'cn'], 'fc');
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs{1}, smmvcnFCs{i}, roiNames, 'cn', ['smmv' num2str(i) 'cn'], 'fc');
%        [adsmadFCsUt, adsmadFCsUtP, adsmadFCsUtP2] = calculateAlzWilcoxonTest(adFCs{1}, smadFCs{i}, roiNames, 'ad', ['sm' num2str(i) 'ad'], 'fc');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs{1}, smcnDLs{i}, roiNames, 'cn', ['sm' num2str(i) 'cn'], 'dlcm');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs{1}, sm2cnDLs{i}, roiNames, 'cn', ['smlin' num2str(i) 'cn'], 'dlcm');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs{1}, smmvcnDLs{i}, roiNames, 'cn', ['smmv' num2str(i) 'cn'], 'dlcm');
%        [adsmadDLsUt, adsmadDLsUtP, adsmadDLsUtP2] = calculateAlzWilcoxonTest(adDLs{1}, smadDLs{i}, roiNames, 'ad', ['sm' num2str(i) 'ad'], 'dlcm');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs{1}, smcnDLWs{i}, roiNames, 'cn', ['sm' num2str(i) 'cn'], 'dlw');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs{1}, sm2cnDLWs{i}, roiNames, 'cn', ['smlin' num2str(i) 'cn'], 'dlw');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs{1}, smmvcnDLWs{i}, roiNames, 'cn', ['smmv' num2str(i) 'cn'], 'dlw');
%        [adsmadDLWsUt, adsmadDLWsUtP, adsmadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs{1}, smadDLWs{i}, roiNames, 'ad', ['sm' num2str(i) 'ad'], 'dlw');
    end
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

    if lags>1, lagStr=num2str(lags); else lagStr=''; end
    if isAutoExo>0, exoStr='_ex'; else exoStr=''; end
    if isempty(activateFunc), linStr='_lin'; else linStr=''; end
    outfName = ['results/adsim3-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

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
            % simulate mVAR network with 1st frame & exogenous input signal
            netMVAR = initMvarNetwork(si, exSignal, [], f.exControl, lags);
            [Y, time] = simulateMvarNetwork(si, exSignal, [], f.exControl, netMVAR);
            [ec, subECs(:,:,i)] = calcDlcmEC(f.netDLCM, [], f.exControl); % not used this
        case 'mpcvarec'
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
            % simulate mPCVAR network with 1st frame & exogenous input signal
            netMVAR = initMpcvarNetwork(si, exSignal, [], f.exControl, lags);
            [Y, time] = simulateMpcvarNetwork(si, exSignal, [], f.exControl, netMVAR);
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
