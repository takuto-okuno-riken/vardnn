
function simulateGenderDLCM
    % CONN fmri data base path :
    base = '../../1000func/';

    % CONN output path
    pathesTR25 = {'nki23-35/2500', 'nki36-45/2500'};
    pathesTR14 = {'nki23-35/1400', 'nki36-45/1400'};
    pathesTR6 = {'nki23-35/645', 'nki36-45/645'};

    % load each type signals
    [tr25.signals, roiNames] = connData2signalsFile(base, pathesTR25, 'tr25', 'data/indi', 'id');
    [tr14.signals] = connData2signalsFile(base, pathesTR14, 'tr14', 'data/indi', 'id');
    [tr6.signals] = connData2signalsFile(base, pathesTR6, 'tr6', 'data/indi', 'id');
    load('data/indi/sbjInfo');
    Idx1 = intersect(find(sbjInfo(:,3)==1),find(sbjInfo(:,4)==1));
    Idx2 = intersect(find(sbjInfo(:,3)==2),find(sbjInfo(:,4)==1));
    Idx3 = intersect(find(sbjInfo2(:,3)==1),find(sbjInfo2(:,4)==1));
    Idx4 = intersect(find(sbjInfo2(:,3)==2),find(sbjInfo2(:,4)==1));

    global resultsPath;
    global resultsPrefix;
    resultsPath = 'results/indi';
    resultsPrefix = 'id';

    maxLag = 5;

    % get DLCM-EC, DLCM-GC & FC from original CN and AD (normal and recovery training)
    tr25 = getDLCM_EC_GC_FC(tr25, roiNames, 'tr25');
    tr14 = getDLCM_EC_GC_FC(tr14, roiNames, 'tr14');
    tr6 = getDLCM_EC_GC_FC(tr6, roiNames, 'tr6');

    % --------------------------------------------------------------------------------------------------------------
    % simulate by algorithms
    tr25 = simulateGroupByAlgorithms(tr25, 'tr25', maxLag);
    tr14 = simulateGroupByAlgorithms(tr14, 'tr14', maxLag);
    tr6  = simulateGroupByAlgorithms(tr6, 'tr6', maxLag);

    % --------------------------------------------------------------------------------------------------------------
    % check DLCM-EC and DLCM-GC of simulated TR25, TR14, TR6
    tr25 = checkSimSignalsByDLCM_EC_GC(tr25, 'tr25', maxLag);
    tr14 = checkSimSignalsByDLCM_EC_GC(tr14, 'tr14', maxLag);
    tr6  = checkSimSignalsByDLCM_EC_GC(tr6, 'tr6', maxLag);

    % check FC of simulated TR25, TR14, TR6
    tr25 = checkSimSignalsByFC(tr25, 'tr25', maxLag);
    tr14 = checkSimSignalsByFC(tr14, 'tr14', maxLag);
    tr6  = checkSimSignalsByFC(tr6, 'tr6', maxLag);

    % --------------------------------------------------------------------------------------------------------------
    % box plot of signal MAEs 
    tr25 = plotSimSignalsResults(tr25, 'tr25', maxLag);

end

% ==================================================================================================================

function g = getDLCM_EC_GC_FC(g, roiNames, groupName)
    signals = g.signals;
    g.sbjNum = length(signals);
    for j=1:1
        % DLCM(i)-GC auto exogenous
        [g.DLs{j}, g.meanDL{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlcm', 0, j, 1);
        % FC no exogenous (pairwise, then exogenous does not have meaning)
        [g.FCs{j}, g.meanFC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'fc', 1);
    end
end

function g = simulateGroupByAlgorithms(g, groupName, maxLag)
    signals = g.signals;
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [g.DLWs{j}, g.smSignals{j}, g.SubDLWs{j}] = simulateNodeSignals(signals, roiNames, groupName, 'dlw', groupName, 0, 0, j, 0);
        g.sigDLWs{j} = (g.DLWs{j} - nanmean(g.DLWs{j}(:))) / nanstd(g.DLWs{j}(:),1);
        g.meanDLW{j} = nanmean(g.DLWs{j},3);

        % DLCM(j) linear no exogenous 
        [g.DLW2s{j}, g.sm2Signals{j}, g.SubDLW2s{j}] = simulateNodeSignals(signals, roiNames, groupName, 'dlw', groupName, 0, 0, j, 0, []);
        g.sigDLW2s{j} = (g.DLW2s{j} - nanmean(g.DLW2s{j}(:))) / nanstd(g.DLW2s{j}(:),1);
        g.meanDLW2{j} = nanmean(g.DLW2s{j},3);

        % mvar(j) no exogenous 
        [g.MVARECs{j}, g.smmvSignals{j}, g.SubMVARECs{j}] = simulateNodeSignals(signals, roiNames, groupName, 'mvarec', groupName, 0, 0, j, 0);
        g.sigMVARECs{j} = (g.MVARECs{j} - nanmean(g.MVARECs{j}(:))) / nanstd(g.MVARECs{j}(:),1);
        g.meanMVAREC{j} = nanmean(g.MVARECs{j},3);

        % mpcvar(j) no exogenous 
        [g.MPCVARECs{j}, g.smmpvSignals{j}, g.SubMPCVARECs{j}] = simulateNodeSignals(signals, roiNames, groupName, 'mpcvarec', groupName, 0, 0, j, 0);
        g.sigMPCVARECs{j} = (g.MPCVARECs{j} - nanmean(g.MPCVARECs{j}(:))) / nanstd(g.MPCVARECs{j}(:),1);
        g.meanMPCVAREC{j} = nanmean(g.MPCVARECs{j},3);
    end

    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [g.DLWs{j}, g.smSignals{j}, g.SubDLWs{j}] = simulateNodeSignals(signals, roiNames, groupName, 'dlw', groupName, 0, 0, i, 1);
        g.sigDLWs{j} = (g.DLWs{j} - nanmean(g.DLWs{j}(:))) / nanstd(g.DLWs{j}(:),1);
        g.meanDLW{j} = nanmean(g.DLWs{j},3);

        % DLCM(j) linear auto exogenous 
        [g.DLW2s{j}, g.sm2Signals{j}, g.SubDLW2s{j}] = simulateNodeSignals(signals, roiNames, groupName, 'dlw', groupName, 0, 0, i, 1, []);
        g.sigDLW2s{j} = (g.DLW2s{j} - nanmean(g.DLW2s{j}(:))) / nanstd(g.DLW2s{j}(:),1);
        g.meanDLW2{j} = nanmean(g.DLW2s{j},3);

        % mvar(j) auto exogenous 
        [g.MVARECs{j}, g.smmvSignals{j}, g.SubMVARECs{j}] = simulateNodeSignals(signals, roiNames, groupName, 'mvarec', groupName, 0, 0, i, 1);
        sigDLWs{j} = (g.MVARECs{j} - nanmean(g.MVARECs{j}(:))) / nanstd(g.MVARECs{j}(:),1);
        meanDLW{j} = nanmean(g.MVARECs{j},3);

        % mpcvar(j) auto exogenous 
        [g.MPCVARECs{j}, g.smmpvSignals{j}, g.SubMPCVARECs{j}] = simulateNodeSignals(signals, roiNames, groupName, 'mpcvarec', groupName, 0, 0, i, 1);
        g.sigMPCVARECs{j} = (g.MPCVARECs{j} - nanmean(g.MPCVARECs{j}(:))) / nanstd(g.MPCVARECs{j}(:),1);
        g.meanMPCVAREC{j} = nanmean(g.MPCVARECs{j},3);
    end
end

function g = checkSimSignalsByDLCM_EC_GC(g, groupName, maxLag)
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [g.smDLs{j}, g.meanSmDL{j}, ~] = calculateConnectivity(g.smSignals{j}, roiNames, ['sm' groupName], 'dlcm', 1, j, 0);
        [g.smDLWs{j}, g.meanSmDLW{j}, ~, g.smSubDLWs{j}] = calculateConnectivity(g.smSignals{j}, roiNames, ['sm' groupName], 'dlw', 1, j, 0);
        % DLCM(j) linear no exogenous 
        [g.sm2DLs{j}, g.meanSm2DL{j}, ~] = calculateConnectivity(g.sm2Signals{j}, roiNames, ['sm' groupName], 'dlcm', 1, j, 0, []);
        [g.sm2DLWs{j}, g.meanSm2DLW{j}, ~, g.sm2SubDLWs{j}] = calculateConnectivity(g.sm2Signals{j}, roiNames, ['sm' groupName], 'dlw', 1, j, 0, []);
        % mvar(j) no exogenous 
        [g.smmvDLs{j}, g.meanSmmvDL{j}, ~] = calculateConnectivity(g.smmvSignals{j}, roiNames, ['smmv' groupName], 'dlcm', 1, j, 0);
        [g.smmvDLWs{j}, g.meanSmmvDLW{j}, ~, g.smmvSubDLW2s{j}] = calculateConnectivity(g.smmvSignals{j}, roiNames, ['smmv' groupName], 'dlw', 1, j, 0);
        % mvar(j) no exogenous 
        [g.smmpvDLs{j}, g.meanSmmpvDL{j}, ~] = calculateConnectivity(g.smmpvSignals{j}, roiNames, ['smmpv' groupName], 'dlcm', 1, j, 0);
        [g.smmpvDLWs{j}, g.meanSmmpvDLW{j}, ~, g.smmpvSubDLW2s{j}] = calculateConnectivity(g.smmpvSignals{j}, roiNames, ['smmpv' groupName], 'dlw', 1, j, 0);
    end
    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [g.smDLs{j}, g.meanSmDL{j}, ~] = calculateConnectivity(g.smSignals{j}, roiNames, ['sm' groupName], 'dlcm', 1, i, 1);
        [g.smDLWs{j}, g.meanSmDLW{j}, ~, g.smSubDLWs{j}] = calculateConnectivity(g.smSignals{j}, roiNames, ['sm' groupName], 'dlw', 1, i, 1);
        % DLCM(j) linear auto exogenous 
        [g.sm2DLs{j}, g.meanSm2DL{j}, ~] = calculateConnectivity(g.sm2Signals{j}, roiNames, ['sm' groupName], 'dlcm', 1, i, 1, []);
        [g.sm2DLWs{j}, g.meanSm2DLW{j}, ~, g.sm2SubDLWs{j}] = calculateConnectivity(g.sm2Signals{j}, roiNames, ['sm' groupName], 'dlw', 1, i, 1, []);
        % mvar(j) auto exogenous 
        [g.smmvDLs{j}, g.meanSmmvDL{j}, ~] = calculateConnectivity(g.smmvSignals{j}, roiNames, ['smmv' groupName], 'dlcm', 1, i, 1);
        [g.smmvDLWs{j}, g.meanSmmvDLW{j}, ~, g.smSubDLW2s{j}] = calculateConnectivity(g.smmvSignals{j}, roiNames, ['smmv' groupName], 'dlw', 1, i, 1);
        % mvar(j) auto exogenous 
        [g.smmpvDLs{j}, g.meanSmmpvDL{j}, ~] = calculateConnectivity(g.smmpvSignals{j}, roiNames, ['smmpv' groupName], 'dlcm', 1, i, 1);
        [g.smmpvDLWs{j}, g.meanSmmpvDLW{j}, ~, g.smmpvSubDLW2s{j}] = calculateConnectivity(g.smmpvSignals{j}, roiNames, ['smmpv' groupName], 'dlw', 1, i, 1);
    end
end

function g = checkSimSignalsByFC(g, groupName, maxLag)
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [g.smFCs{j}, g.meanSmFC{j}, ~] = calculateConnectivity(g.smSignals{j}, roiNames, ['sm' groupName], 'fc', 1, j, 0);
        % DLCM(j) no exogenous 
        [g.sm2FCs{j}, g.meanSm2FC{j}, ~] = calculateConnectivity(g.sm2Signals{j}, roiNames, ['sm' groupName], 'fc', 1, j, 0, []);
        % mvar(j) no exogenous 
        [g.smmvFCs{j}, g.meanSmmvFC{j}, ~] = calculateConnectivity(g.smmvSignals{j}, roiNames, ['smmv' groupName], 'fc', 1, j, 0);
        % mpcvar(j) no exogenous 
        [g.smmpvFCs{j}, g.meanSmmpvFC{j}, ~] = calculateConnectivity(g.smmpvSignals{j}, roiNames, ['smmpv' groupName], 'fc', 1, j, 0);
    end
    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [g.smFCs{j}, g.meanSmFC{j}, ~] = calculateConnectivity(g.smSignals{j}, roiNames, ['sm' groupName], 'fc', 1, i, 1);
        % DLCM(j) auto exogenous 
        [g.sm2FCs{j}, g.meanSm2FC{j}, ~] = calculateConnectivity(g.sm2Signals{j}, roiNames, ['sm' groupName], 'fc', 1, i, 1, []);
        % mvar(j) auto exogenous 
        [g.smmvFCs{j}, g.meanSmmvFC{j}, ~] = calculateConnectivity(g.smmvSignals{j}, roiNames, ['smmv' groupName], 'fc', 1, i, 1);
        % mpcvar(j) auto exogenous 
        [g.smmpvFCs{j}, g.meanSmmpvFC{j}, ~] = calculateConnectivity(g.smmpvSignals{j}, roiNames, ['smmpv' groupName], 'fc', 1, i, 1);
    end
end

function g = plotSimSignalsResults(g, groupName, maxLag)
    signals = g.signals;
    maes = nan(g.sbjNum,40);
    for j=1:g.sbjNum
        si = convert2SigmoidSignal(signals{j});
        for i=1:maxLag*2
            maes(j,i) = getTwoSignalsError(si, g.smSignals{i}{j}); % non linear (no ex, ex)
            maes(j,10+i) = getTwoSignalsError(si, g.sm2Signals{i}{j}); % linear (no ex, ex)
            maes(j,20+i) = getTwoSignalsError(si, g.smmvSignals{i}{j}); % linear (no ex, ex)
            maes(j,30+i) = getTwoSignalsError(si, g.smmpvSignals{i}{j}); % linear (no ex, ex)
        end
    end
    figure; boxplot(maes); title(['MAEs between mean cnSignals and each algorithm signals : ' groupName]);

    % plot correlation and cos similarity
    cosSim = zeros(120,1);
    for k=1:maxLag*2
        i=k;
        cosSim(i) = getCosSimilarity(g.meanDLW{1}, g.meanSmDLW{k}); i=i+10; % non linear (no ex, ex)
        cosSim(i) = getCosSimilarity(g.meanDLW{1}, g.meanSm2DLW{k}); i=i+10; % linear (no ex, ex)
        cosSim(i) = getCosSimilarity(g.meanDLW{1}, g.meanSmmvDLW{k}); i=i+10; % mvar linear (no ex, ex)
        cosSim(i) = getCosSimilarity(g.meanDLW{1}, g.meanSmmpvDLW{k}); i=i+10; % mpcvar linear (no ex, ex)
        cosSim(i) = getCosSimilarity(g.meanDL{1}, g.meanSmDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanDL{1}, g.meanSm2DL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanDL{1}, g.meanSmmvDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanDL{1}, g.meanSmmpvDL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanFC{1}, g.meanSmFC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanFC{1}, g.meanSm2FC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanFC{1}, g.meanSmmvFC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(g.meanFC{1}, g.meanSmmpvFC{k}); i=i+10;
    end
    figure; bar(cosSim); title(['cos similarity between mean CN matrix and SimCN by each algorithm : ' groupName]);
%{
    cosSim = zeros(algNum,1);
    figure; bar(cosSim); title('cos similarity between mean AD matrix and SimAD by each algorithm');
%}  
    % plot box-and-whisker plot of cos similarity between mean ec matrix and each subject ec
    cosSims = nan(g.sbjNum,120);
    for j=1:g.sbjNum
        for k=1:maxLag*2
            i=k;
            cosSims(j,i) = getCosSimilarity(g.meanDLW{1}, g.smDLWs{k}(:,:,j)); i=i+10; % non linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.meanDLW{1}, g.sm2DLWs{k}(:,:,j)); i=i+10; % linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.meanDLW{1}, g.smmvDLWs{k}(:,:,j)); i=i+10; % mvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.meanDLW{1}, g.smmpvDLWs{k}(:,:,j)); i=i+10; % mpcvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.meanDL{1}, g.smDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanDL{1}, g.sm2DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanDL{1}, g.smmvDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanDL{1}, g.smmpvDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanFC{1}, g.smFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanFC{1}, g.sm2FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanFC{1}, g.smmvFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.meanFC{1}, g.smmpvFCs{k}(:,:,j)); i=i+10;
        end
%        cosSims(j,1) = getCosSimilarity(meanCnDLW, cnDLWs(:,:,i));
%        cosSims(i,11) = getCosSimilarity(meanCnDL, cnDLs(:,:,i));
%        cosSims(i,21) = getCosSimilarity(meanCnFC, cnFCs(:,:,i));
    end
    figure; boxplot(cosSims); title(['cos similarity between mean CN matrix and each subject EC by each algorithm : ' groupName]);
%    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:30));
%    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
%    [sortCosCnDLW,idxCosCnDLW] = sort(cosSims(:,1),'descend');
%    [sortCosCnFC,idxCosCnFC] = sort(cosSims(:,21),'descend');

    % plot correlation between original EC matrix and simulated signals EC matrix
%{
    for i=1:1 %cnSbjNum
        % calc & plot correlation of original Zi vs simulating Zi
        [corrZi, corrZij] = calcCorrelationZiZij(cnSubDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum);
        plotCorrelationZiZij(cnDLWs(:,:,i), cnSubDLWs(:,:,i), smcnDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum, ['sbj' num2str(i)], 'original', 'simulating');
    end
%}
    % plot box-and-whisker plot of cos similarity between original ec matrix and simulated signals ec matrix
    cosSims = nan(g.sbjNum,120);
    for j=1:g.sbjNum
        for k=1:maxLag*2
            i=k;
            cosSims(j,i) = getCosSimilarity(g.DLWs{1}(:,:,j), g.smDLWs{k}(:,:,j)); i=i+10; % non linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.DLWs{1}(:,:,j), g.sm2DLWs{k}(:,:,j)); i=i+10; % linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.DLWs{1}(:,:,j), g.smmvDLWs{k}(:,:,j)); i=i+10; % mvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.DLWs{1}(:,:,j), g.smmpvDLWs{k}(:,:,j)); i=i+10; % mpcvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(g.DLs{1}(:,:,j), g.smDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.DLs{1}(:,:,j), g.sm2DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.DLs{1}(:,:,j), g.smmvDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.DLs{1}(:,:,j), g.smmpvDLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.FCs{1}(:,:,j), g.smFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.FCs{1}(:,:,j), g.sm2FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.FCs{1}(:,:,j), g.smmvFCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(g.FCs{1}(:,:,j), g.smmpvFCs{k}(:,:,j)); i=i+10;
        end
    end
    figure; boxplot(cosSims); title(['cos similarity between mean CN matrix and each subject EC by each algorithm : ' groupName]);

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
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(g.FCs{1}, g.smFCs{i}, roiNames, groupName, ['sm' num2str(i) groupName], 'fc');
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(g.FCs{1}, g.sm2FCs{i}, roiNames, groupName, ['smlin' num2str(i) groupName], 'fc');
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(g.FCs{1}, g.smmvFCs{i}, roiNames, groupName, ['smmv' num2str(i) groupName], 'fc');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(g.DLs{1}, g.smDLs{i}, roiNames, groupName, ['sm' num2str(i) groupName], 'dlcm');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(g.DLs{1}, g.sm2DLs{i}, roiNames, groupName, ['smlin' num2str(i) groupName], 'dlcm');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(g.DLs{1}, g.smmvDLs{i}, roiNames, groupName, ['smmv' num2str(i) groupName], 'dlcm');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(g.DLWs{1}, g.smDLWs{i}, roiNames, groupName, ['sm' num2str(i) groupName], 'dlw');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(g.DLWs{1}, g.sm2DLWs{i}, roiNames, groupName, ['smlin' num2str(i) groupName], 'dlw');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(g.DLWs{1}, g.smmvDLWs{i}, roiNames, groupName, ['smmv' num2str(i) groupName], 'dlw');
    end
end

% ---------------------------------------------------------------------------------------------------------------------------------------

function [ECs, simSignals, subECs] = simulateNodeSignals(signals, roiNames, group, algorithm, orgGroup, isRaw, inSiRange, lags, isAutoExo, activateFunc)
    if nargin < 10, activateFunc = @reluLayer; end
    if nargin < 9, isAutoExo = 0; end
    if nargin < 8, lags = 1; end
    if nargin < 7, inSiRange = 0; end
    if nargin < 6, isRaw = 0; end

    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    global resultsPath;
    global resultsPrefix;
    if lags>1, lagStr=num2str(lags); else lagStr=''; end
    if isAutoExo>0, exoStr='_ex'; else exoStr=''; end
    if isempty(activateFunc), linStr='_lin'; else linStr=''; end
    outfName = [resultsPath 'sim/' resultsPrefix 'sim-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '.mat'];
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
            dlcmName = [resultsPath '/' resultsPrefix '-dlcm' lagStr exoStr linStr '-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
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
            dlcmName = [resultsPath '/' resultsPrefix '-dlcm' lagStr exoStr linStr '-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
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
            dlcmName = [resultsPath '/' resultsPrefix '-dlcm' lagStr exoStr linStr '-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
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

    global resultsPath;
    global resultsPrefix;
    outfName = [resultsPath 'sim/' resultsPrefix 'sim-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROWNUM) '-utest.mat'];
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
