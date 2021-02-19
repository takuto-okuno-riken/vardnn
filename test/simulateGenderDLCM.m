
function simulateGenderDLCM
    % CONN fmri data base path :
    base = '../../1000func/';

    % CONN output path
    pathesTR25 = {'nki23-35/2500', 'nki36-45/2500'};
    pathesTR14 = {'nki23-35/1400', 'nki36-45/1400'};
    pathesTR6 = {'nki23-35/645', 'nki36-45/645'};

    % load each type signals
    [tr25Signals, roiNames] = connData2signalsFile(base, pathesTR25, 'tr25', 'data/indi', 'id');
    [tr14Signals] = connData2signalsFile(base, pathesTR14, 'tr14', 'data/indi', 'id');
    [tr6Signals] = connData2signalsFile(base, pathesTR6, 'tr6', 'data/indi', 'id');
    load('data/indi/sbjInfo');
    Idx1 = intersect(find(sbjInfo(:,3)==1),find(sbjInfo(:,4)==1));
    Idx2 = intersect(find(sbjInfo(:,3)==2),find(sbjInfo(:,4)==1));
    Idx3 = intersect(find(sbjInfo2(:,3)==1),find(sbjInfo2(:,4)==1));
    Idx4 = intersect(find(sbjInfo2(:,3)==2),find(sbjInfo2(:,4)==1));

    global resultsPath;
    global resultsPrefix;
    resultsPath = 'results/indi';
    resultsPrefix = 'id';

    tr25SbjNum = length(tr25Signals);
    tr14SbjNum = length(tr14Signals);
    nodeNum = size(tr25Signals{1},1);
    maxLag = 5;

    % get DLCM-EC, DLCM-GC & FC from original CN and AD (normal and recovery training)
    for j=1:1
        % DLCM(i)-GC auto exogenous
        [tr25DLs{j}, meanTr25DL{j}, ~] = calculateConnectivity(tr25Signals, roiNames, 'tr25', 'dlcm', 0, j, 1);
        [tr14DLs{j}, meanTr14DL{j}, ~] = calculateConnectivity(tr14Signals, roiNames, 'tr14', 'dlcm', 0, j, 1);
        [tr6DLs{j}, meanTr6DL{j}, ~] = calculateConnectivity(tr6Signals, roiNames, 'tr6', 'dlcm', 0, j, 1);
        % FC no exogenous (pairwise, then exogenous does not have meaning)
        [tr25FCs{j}, meanTr25FC{j}, ~] = calculateConnectivity(tr25Signals, roiNames, 'tr25', 'fc', 1);
        [tr14FCs{j}, meanTr14FC{j}, ~] = calculateConnectivity(tr14Signals, roiNames, 'tr14', 'fc', 1);
        [tr6FCs{j}, meanTr6FC{j}, ~] = calculateConnectivity(tr6Signals, roiNames, 'tr6', 'fc', 1);
    end

    % --------------------------------------------------------------------------------------------------------------
    % simulate by algorithms
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [tr25DLWs{j}, smtr25Signals{j}, tr25SubDLWs{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'dlw', 'tr25', 0, 0, j, 0);
        [tr14DLWs{j}, smtr14Signals{j}, tr14SubDLWs{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'dlw', 'tr14', 0, 0, j, 0);
        sigTr25DLWs{j} = (tr25DLWs{j} - nanmean(tr25DLWs{j}(:))) / nanstd(tr25DLWs{j}(:),1);
        sigTr14DLWs{j} = (tr14DLWs{j} - nanmean(tr14DLWs{j}(:))) / nanstd(tr14DLWs{j}(:),1);
        meanTr25DLW{j} = nanmean(tr25DLWs{j},3);
        meanTr14DLW{j} = nanmean(tr14DLWs{j},3);

        % DLCM(j) linear no exogenous 
        [tr25DLW2s{j}, sm2tr25Signals{j}, tr25SubDLW2s{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'dlw', 'tr25', 0, 0, j, 0, []);
        [tr14DLW2s{j}, sm2tr14Signals{j}, tr14SubDLW2s{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'dlw', 'tr14', 0, 0, j, 0, []);
        sigTr25DLW2s{j} = (tr25DLW2s{j} - nanmean(tr25DLW2s{j}(:))) / nanstd(tr25DLW2s{j}(:),1);
        sigTr14DLW2s{j} = (tr14DLW2s{j} - nanmean(tr14DLW2s{j}(:))) / nanstd(tr14DLW2s{j}(:),1);
        meanTr25DLW2{j} = nanmean(tr25DLW2s{j},3);
        meanTr14DLW2{j} = nanmean(tr14DLW2s{j},3);

        % mvar(j) no exogenous 
        [tr25MVARECs{j}, smmvtr25Signals{j}, tr25SubMVARECs{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'mvarec', 'tr25', 0, 0, j, 0);
        [tr14MVARECs{j}, smmvtr14Signals{j}, tr14SubMVARECs{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'mvarec', 'tr14', 0, 0, j, 0);
        sigTr25MVARECs{j} = (tr25MVARECs{j} - nanmean(tr25MVARECs{j}(:))) / nanstd(tr25MVARECs{j}(:),1);
        sigTr14MVARECs{j} = (tr14MVARECs{j} - nanmean(tr14MVARECs{j}(:))) / nanstd(tr14MVARECs{j}(:),1);
        meanTr25MVAREC{j} = nanmean(tr25MVARECs{j},3);
        meanTr14MVAREC{j} = nanmean(tr14MVARECs{j},3);

        % mpcvar(j) no exogenous 
        [tr25MPCVARECs{j}, smmpvtr25Signals{j}, tr25SubMPCVARECs{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'mpcvarec', 'tr25', 0, 0, j, 0);
        [tr14MPCVARECs{j}, smmpvtr25Signals{j}, tr14SubMPCVARECs{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'mpcvarec', 'tr14', 0, 0, j, 0);
        sigTr25MPCVARECs{j} = (tr25MPCVARECs{j} - nanmean(tr25MPCVARECs{j}(:))) / nanstd(tr25MPCVARECs{j}(:),1);
        sigTr14MPCVARECs{j} = (tr14MPCVARECs{j} - nanmean(tr14MPCVARECs{j}(:))) / nanstd(tr14MPCVARECs{j}(:),1);
        meanTr25MPCVAREC{j} = nanmean(tr25MPCVARECs{j},3);
        meanTr14MPCVAREC{j} = nanmean(tr14MPCVARECs{j},3);
    end

    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [tr25DLWs{j}, smtr25Signals{j}, tr25SubDLWs{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'dlw', 'tr25', 0, 0, i, 1);
        [tr14DLWs{j}, smtr14Signals{j}, tr14SubDLWs{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'dlw', 'tr14', 0, 0, i, 1);
        sigTr25DLWs{j} = (tr25DLWs{j} - nanmean(tr25DLWs{j}(:))) / nanstd(tr25DLWs{j}(:),1);
        sigTr14DLWs{j} = (tr14DLWs{j} - nanmean(tr14DLWs{j}(:))) / nanstd(tr14DLWs{j}(:),1);
        meanTr25DLW{j} = nanmean(tr25DLWs{j},3);
        meanTr14DLW{j} = nanmean(tr14DLWs{j},3);

        % DLCM(j) linear auto exogenous 
        [tr25DLW2s{j}, sm2tr25Signals{j}, tr25SubDLW2s{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'dlw', 'tr25', 0, 0, i, 1, []);
        [tr14DLW2s{j}, sm2tr14Signals{j}, tr14SubDLW2s{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'dlw', 'tr14', 0, 0, i, 1, []);
        sigTr25DLW2s{j} = (tr25DLW2s{j} - nanmean(tr25DLW2s{j}(:))) / nanstd(tr25DLW2s{j}(:),1);
        sigTr14DLW2s{j} = (tr14DLW2s{j} - nanmean(tr14DLW2s{j}(:))) / nanstd(tr14DLW2s{j}(:),1);
        meanTr25DLW2{j} = nanmean(tr25DLW2s{j},3);
        meanTr14DLW2{j} = nanmean(tr14DLW2s{j},3);

        % mvar(j) auto exogenous 
        [tr25MVARECs{j}, smmvtr25Signals{j}, tr25SubMVARECs{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'mvarec', 'tr25', 0, 0, i, 1);
        [tr14MVARECs{j}, smmvtr14Signals{j}, tr14SubMVARECs{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'mvarec', 'tr14', 0, 0, i, 1);
        sigTr25DLWs{j} = (tr25MVARECs{j} - nanmean(tr25MVARECs{j}(:))) / nanstd(tr25MVARECs{j}(:),1);
        sigTr14DLWs{j} = (tr14MVARECs{j} - nanmean(tr14MVARECs{j}(:))) / nanstd(tr14MVARECs{j}(:),1);
        meanTr25DLW{j} = nanmean(tr25MVARECs{j},3);
        meanTr14DLW{j} = nanmean(tr14MVARECs{j},3);

        % mpcvar(j) auto exogenous 
        [tr25MPCVARECs{j}, smmpvtr25Signals{j}, tr25SubMPCVARECs{j}] = simulateNodeSignals(tr25Signals, roiNames, 'tr25', 'mpcvarec', 'tr25', 0, 0, i, 1);
        [tr14MPCVARECs{j}, smmpvtr14Signals{j}, tr14SubMPCVARECs{j}] = simulateNodeSignals(tr14Signals, roiNames, 'tr14', 'mpcvarec', 'tr14', 0, 0, i, 1);
        sigTr25MPCVARECs{j} = (tr25MPCVARECs{j} - nanmean(tr25MPCVARECs{j}(:))) / nanstd(tr25MPCVARECs{j}(:),1);
        sigTr14MPCVARECs{j} = (tr14MPCVARECs{j} - nanmean(tr14MPCVARECs{j}(:))) / nanstd(tr14MPCVARECs{j}(:),1);
        meanTr25MPCVAREC{j} = nanmean(tr25MPCVARECs{j},3);
        meanTr14MPCVAREC{j} = nanmean(tr14MPCVARECs{j},3);
    end

    % --------------------------------------------------------------------------------------------------------------
    % check DLCM-EC and DLCM-GC of simulated CN and AD
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [smtr25DLs{j}, meanSmtr25DL{j}, ~] = calculateConnectivity(smtr25Signals{j}, roiNames, 'smtr25', 'dlcm', 1, j, 0);
        [smtr25DLWs{j}, meanSmtr25DLW{j}, ~, smtr25SubDLWs{j}] = calculateConnectivity(smtr25Signals{j}, roiNames, 'smtr25', 'dlw', 1, j, 0);
        % DLCM(j) linear no exogenous 
        [sm2tr25DLs{j}, meanSm2tr25DL{j}, ~] = calculateConnectivity(sm2tr25Signals{j}, roiNames, 'smtr25', 'dlcm', 1, j, 0, []);
        [sm2tr25DLWs{j}, meanSm2tr25DLW{j}, ~, sm2tr25SubDLWs{j}] = calculateConnectivity(sm2tr25Signals{j}, roiNames, 'smtr25', 'dlw', 1, j, 0, []);
        % mvar(j) no exogenous 
        [smmvtr25DLs{j}, meanSmmvtr25DL{j}, ~] = calculateConnectivity(smmvtr25Signals{j}, roiNames, 'smmvtr25', 'dlcm', 1, j, 0);
        [smmvtr25DLWs{j}, meanSmmvtr25DLW{j}, ~, smmvtr25SubDLW2s{j}] = calculateConnectivity(smmvtr25Signals{j}, roiNames, 'smmvtr25', 'dlw', 1, j, 0);
        % mvar(j) no exogenous 
        [smmpvtr25DLs{j}, meanSmmpvtr25DL{j}, ~] = calculateConnectivity(smmpvtr25Signals{j}, roiNames, 'smmpvtr25', 'dlcm', 1, j, 0);
        [smmpvtr25DLWs{j}, meanSmmpvtr25DLW{j}, ~, smmpvtr25SubDLW2s{j}] = calculateConnectivity(smmpvtr25Signals{j}, roiNames, 'smmpvtr25', 'dlw', 1, j, 0);
    end
    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [smtr25DLs{j}, meanSmtr25DL{j}, ~] = calculateConnectivity(smtr25Signals{j}, roiNames, 'smtr25', 'dlcm', 1, i, 1);
        [smtr25DLWs{j}, meanSmtr25DLW{j}, ~, smtr25SubDLWs{j}] = calculateConnectivity(smtr25Signals{j}, roiNames, 'smtr25', 'dlw', 1, i, 1);
        % DLCM(j) linear auto exogenous 
        [sm2tr25DLs{j}, meanSm2tr25DL{j}, ~] = calculateConnectivity(sm2tr25Signals{j}, roiNames, 'smtr25', 'dlcm', 1, i, 1, []);
        [sm2tr25DLWs{j}, meanSm2tr25DLW{j}, ~, sm2tr25SubDLWs{j}] = calculateConnectivity(sm2tr25Signals{j}, roiNames, 'smtr25', 'dlw', 1, i, 1, []);
        % mvar(j) auto exogenous 
        [smmvtr25DLs{j}, meanSmmvtr25DL{j}, ~] = calculateConnectivity(smmvtr25Signals{j}, roiNames, 'smmvtr25', 'dlcm', 1, i, 1);
        [smmvtr25DLWs{j}, meanSmmvtr25DLW{j}, ~, smtr25SubDLW2s{j}] = calculateConnectivity(smmvtr25Signals{j}, roiNames, 'smmvtr25', 'dlw', 1, i, 1);
        % mvar(j) auto exogenous 
        [smmpvtr25DLs{j}, meanSmmpvtr25DL{j}, ~] = calculateConnectivity(smmpvtr25Signals{j}, roiNames, 'smmpvtr25', 'dlcm', 1, i, 1);
        [smmpvtr25DLWs{j}, meanSmmpvtr25DLW{j}, ~, smmpvtr25SubDLW2s{j}] = calculateConnectivity(smmpvtr25Signals{j}, roiNames, 'smmpvtr25', 'dlw', 1, i, 1);
    end

    % check FC of simulated CN and AD
    for j=1:maxLag
        % DLCM(j) no exogenous 
        [smtr25FCs{j}, meanSmtr25FC{j}, ~] = calculateConnectivity(smtr25Signals{j}, roiNames, 'smtr25', 'fc', 1, j, 0);
        [smtr14FCs{j}, meanSmtr14FC{j}, ~] = calculateConnectivity(smtr14Signals{j}, roiNames, 'smtr14', 'fc', 1, j, 0);
        % DLCM(j) no exogenous 
        [sm2tr25FCs{j}, meanSm2tr25FC{j}, ~] = calculateConnectivity(sm2tr25Signals{j}, roiNames, 'smtr25', 'fc', 1, j, 0, []);
        [sm2tr14FCs{j}, meanSm2tr14FC{j}, ~] = calculateConnectivity(sm2tr14Signals{j}, roiNames, 'smtr14', 'fc', 1, j, 0, []);
        % mvar(j) no exogenous 
        [smmvtr25FCs{j}, meanSmmvtr25FC{j}, ~] = calculateConnectivity(smmvtr25Signals{j}, roiNames, 'smmvtr25', 'fc', 1, j, 0);
        [smmvtr14FCs{j}, meanSmmvtr14FC{j}, ~] = calculateConnectivity(smmvtr14Signals{j}, roiNames, 'smmvtr14', 'fc', 1, j, 0);
        % mpcvar(j) no exogenous 
        [smmpvtr25FCs{j}, meanSmmpvtr25FC{j}, ~] = calculateConnectivity(smmpvtr25Signals{j}, roiNames, 'smmpvtr25', 'fc', 1, j, 0);
        [smmpvtr14FCs{j}, meanSmmpvtr14FC{j}, ~] = calculateConnectivity(smmpvtr14Signals{j}, roiNames, 'smmpvtr14', 'fc', 1, j, 0);
    end
    for i=1:maxLag
        j = i+maxLag;
        % DLCM(j) auto exogenous 
        [smtr25FCs{j}, meanSmtr25FC{j}, ~] = calculateConnectivity(smtr25Signals{j}, roiNames, 'smtr25', 'fc', 1, i, 1);
        [smtr14FCs{j}, meanSmtr14FC{j}, ~] = calculateConnectivity(smtr14Signals{j}, roiNames, 'smtr14', 'fc', 1, i, 1);
        % DLCM(j) auto exogenous 
        [sm2tr25FCs{j}, meanSm2tr25FC{j}, ~] = calculateConnectivity(sm2tr25Signals{j}, roiNames, 'smtr25', 'fc', 1, i, 1, []);
        [sm2tr14FCs{j}, meanSm2tr14FC{j}, ~] = calculateConnectivity(sm2tr14Signals{j}, roiNames, 'smtr14', 'fc', 1, i, 1, []);
        % mvar(j) auto exogenous 
        [smmvtr25FCs{j}, meanSmmvtr25FC{j}, ~] = calculateConnectivity(smmvtr25Signals{j}, roiNames, 'smmvtr25', 'fc', 1, i, 1);
        [smmvtr14FCs{j}, meanSmmvtr14FC{j}, ~] = calculateConnectivity(smmvtr14Signals{j}, roiNames, 'smmvtr14', 'fc', 1, i, 1);
        % mpcvar(j) auto exogenous 
        [smmpvtr25FCs{j}, meanSmmpvtr25FC{j}, ~] = calculateConnectivity(smmpvtr25Signals{j}, roiNames, 'smmpvtr25', 'fc', 1, i, 1);
        [smmpvtr14FCs{j}, meanSmmpvtr14FC{j}, ~] = calculateConnectivity(smmpvtr14Signals{j}, roiNames, 'smmpvtr14', 'fc', 1, i, 1);
    end

    % --------------------------------------------------------------------------------------------------------------
    % box plot of signal MAEs 
    maes = nan(tr25SbjNum,40);
    for j=1:tr25SbjNum
        si = convert2SigmoidSignal(tr25Signals{j});
        for i=1:maxLag*2
            maes(j,i) = getTwoSignalsError(si, smtr25Signals{i}{j}); % non linear (no ex, ex)
            maes(j,10+i) = getTwoSignalsError(si, sm2tr25Signals{i}{j}); % linear (no ex, ex)
            maes(j,20+i) = getTwoSignalsError(si, smmvtr25Signals{i}{j}); % linear (no ex, ex)
            maes(j,30+i) = getTwoSignalsError(si, smmpvtr25Signals{i}{j}); % linear (no ex, ex)
        end
    end
    figure; boxplot(maes); title('MAEs between mean cnSignals and each algorithm signals');

    % plot correlation and cos similarity
    cosSim = zeros(120,1);
    for k=1:maxLag*2
        i=k;
        cosSim(i) = getCosSimilarity(meanTr25DLW{1}, meanSmtr25DLW{k}); i=i+10; % non linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanTr25DLW{1}, meanSm2tr25DLW{k}); i=i+10; % linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanTr25DLW{1}, meanSmmvtr25DLW{k}); i=i+10; % mvar linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanTr25DLW{1}, meanSmmpvtr25DLW{k}); i=i+10; % mpcvar linear (no ex, ex)
        cosSim(i) = getCosSimilarity(meanTr25DL{1}, meanSmtr25DL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25DL{1}, meanSm2tr25DL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25DL{1}, meanSmmvtr25DL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25DL{1}, meanSmmpvtr25DL{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25FC{1}, meanSmtr25FC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25FC{1}, meanSm2tr25FC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25FC{1}, meanSmmvtr25FC{k}); i=i+10;
        cosSim(i) = getCosSimilarity(meanTr25FC{1}, meanSmmpvtr25FC{k}); i=i+10;
    end
    figure; bar(cosSim); title('cos similarity between mean CN matrix and SimCN by each algorithm');
%{
    cosSim = zeros(algNum,1);
    figure; bar(cosSim); title('cos similarity between mean AD matrix and SimAD by each algorithm');
%}  
    % plot box-and-whisker plot of cos similarity between mean ec matrix and each subject ec
    cosSims = nan(tr25SbjNum,120);
    for j=1:tr25SbjNum
        for k=1:maxLag*2
            i=k;
            cosSims(j,i) = getCosSimilarity(meanTr25DLW{1}, smtr25DLWs{k}(:,:,j)); i=i+10; % non linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanTr25DLW{1}, sm2tr25DLWs{k}(:,:,j)); i=i+10; % linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanTr25DLW{1}, smmvtr25DLWs{k}(:,:,j)); i=i+10; % mvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanTr25DLW{1}, smmpvtr25DLWs{k}(:,:,j)); i=i+10; % mpcvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(meanTr25DL{1}, smtr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25DL{1}, sm2tr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25DL{1}, smmvtr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25DL{1}, smmpvtr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25FC{1}, smtr25FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25FC{1}, sm2tr25FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25FC{1}, smmvtr25FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(meanTr25FC{1}, smmpvtr25FCs{k}(:,:,j)); i=i+10;
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

    cosSims = nan(tr14SbjNum,90);
    for j=1:tr14SbjNum
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
    cosSims = nan(tr25SbjNum,120);
    for j=1:tr25SbjNum
        for k=1:maxLag*2
            i=k;
            cosSims(j,i) = getCosSimilarity(tr25DLWs{1}(:,:,j), smtr25DLWs{k}(:,:,j)); i=i+10; % non linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(tr25DLWs{1}(:,:,j), sm2tr25DLWs{k}(:,:,j)); i=i+10; % linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(tr25DLWs{1}(:,:,j), smmvtr25DLWs{k}(:,:,j)); i=i+10; % mvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(tr25DLWs{1}(:,:,j), smmpvtr25DLWs{k}(:,:,j)); i=i+10; % mpcvar linear (no ex, ex)
            cosSims(j,i) = getCosSimilarity(tr25DLs{1}(:,:,j), smtr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25DLs{1}(:,:,j), sm2tr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25DLs{1}(:,:,j), smmvtr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25DLs{1}(:,:,j), smmpvtr25DLs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25FCs{1}(:,:,j), smtr25FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25FCs{1}(:,:,j), sm2tr25FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25FCs{1}(:,:,j), smmvtr25FCs{k}(:,:,j)); i=i+10;
            cosSims(j,i) = getCosSimilarity(tr25FCs{1}(:,:,j), smmpvtr25FCs{k}(:,:,j)); i=i+10;
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
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(tr25FCs{1}, smtr25FCs{i}, roiNames, 'tr25', ['sm' num2str(i) 'tr25'], 'fc');
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(tr25FCs{1}, sm2tr25FCs{i}, roiNames, 'tr25', ['smlin' num2str(i) 'tr25'], 'fc');
        [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(tr25FCs{1}, smmvtr25FCs{i}, roiNames, 'tr25', ['smmv' num2str(i) 'tr25'], 'fc');
%        [adsmadFCsUt, adsmadFCsUtP, adsmadFCsUtP2] = calculateAlzWilcoxonTest(adFCs{1}, smadFCs{i}, roiNames, 'tr14', ['sm' num2str(i) 'tr14'], 'fc');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(tr25DLs{1}, smtr25DLs{i}, roiNames, 'tr25', ['sm' num2str(i) 'tr25'], 'dlcm');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(tr25DLs{1}, sm2tr25DLs{i}, roiNames, 'tr25', ['smlin' num2str(i) 'tr25'], 'dlcm');
        [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(tr25DLs{1}, smmvtr25DLs{i}, roiNames, 'tr25', ['smmv' num2str(i) 'tr25'], 'dlcm');
%        [adsmadDLsUt, adsmadDLsUtP, adsmadDLsUtP2] = calculateAlzWilcoxonTest(adDLs{1}, smadDLs{i}, roiNames, 'tr14', ['sm' num2str(i) 'tr14'], 'dlcm');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(tr25DLWs{1}, smtr25DLWs{i}, roiNames, 'tr25', ['sm' num2str(i) 'tr25'], 'dlw');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(tr25DLWs{1}, sm2tr25DLWs{i}, roiNames, 'tr25', ['smlin' num2str(i) 'tr25'], 'dlw');
        [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(tr25DLWs{1}, smmvtr25DLWs{i}, roiNames, 'tr25', ['smmv' num2str(i) 'tr25'], 'dlw');
%        [adsmadDLWsUt, adsmadDLWsUtP, adsmadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs{1}, smadDLWs{i}, roiNames, 'tr14', ['sm' num2str(i) 'tr14'], 'dlw');
    end
end

% ==================================================================================================================

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
