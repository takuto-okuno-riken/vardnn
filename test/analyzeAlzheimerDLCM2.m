function analyzeAlzheimerDLCM2
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};
%    pathesMCI = {'ADNI2_65-75_F_MCI_nii', 'ADNI2_65-75_M_MCI_nii'};

    % load each type signals
    [cn.signals, roiNames] = connData2signalsFile(base, pathesCN, 'cn', 'data/ad', 'ad');
    [ad.signals] = connData2signalsFile(base, pathesAD, 'ad', 'data/ad', 'ad');
%    [mciSignals] = connData2signalsFile(base, pathesMCI, 'mci', 'data/ad', 'ad');

    global resultsPath;
    global resultsPrefix;
    resultsPath = 'results/ad';
    resultsPrefix = 'ad';

    maxLag = 5;
    
    % calculate connectivity
    cn = calculateConnectivitiesByAlgorithms(cn, roiNames, 'cn', maxLag);
    ad = calculateConnectivitiesByAlgorithms(ad, roiNames, 'ad', maxLag);
    
    % diagnose groups and show ROC curves
    statisticalGroupIdentificationByAlgorithms(cn, ad, roiNames, 'cn', 'ad', maxLag);
end


function g = calculateConnectivitiesByAlgorithms(g, roiNames, groupName, maxLag)
    signals = g.signals;
    for j=1:maxLag
        % mvGC(i) no exogenous 
        [g.GCs{j}, g.meanGC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'gc', 1, j, 0);
        % mvarEC(i) no exogenous 
        [g.MVARECs{j}, g.meanMVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mvarec', 1, j, 0);
        [g.MVARs{j}, g.meanMVAR{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mvar', 1, j, 0);
        % pvarEC(i) no exogenous 
        [g.PVARECs{j}, g.meanPVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'pvarec', 1, j, 0);
        % mpcvarEC(i) no exogenous 
        [g.MPCVARECs{j}, g.meanMPCVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mpcvarec', 1, j, 0);
        % ppcvarEC(i) no exogenous 
        [g.PPCVARECs{j}, g.meanPPCVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'ppcvarec', 1, j, 0);
        % mpcvarGC(i) no exogenous 
        [g.MPCVARGCs{j}, g.meanMPCVARGC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mpcvargc', 1, j, 0);
        % ppcvarGC(i) no exogenous 
        [g.PPCVARGCs{j}, g.meanPPCVARGC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'ppcvargc', 1, j, 0);
        % DLCM(i)-GC linear no exogenous
        [g.DL2s{j}, g.meanDL2{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlcm', 0, j, 0, []);
        % DLCM(i)-EC linear no exogenous
        [g.DLW2s{j}, g.meanDLW2{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlw', 0, j, 0, []);
        % DLCM(i)-GC no exogenous
        [g.DLs{j}, g.meanDL{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlcm', 0, j, 0);
        % DLCM(i)-EC no exogenous
        [g.DLWs{j}, g.meanDLW{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlw', 0, j, 0);
    end

    for i=1:maxLag
        j = i+maxLag;
        % mvGC(i) auto exogenous 
        [g.GCs{j}, g.meanGC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'gc', 1, i, 1);
        % mvarEC(i) auto exogenous 
        [g.MVARECs{j}, g.meanMVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mvarec', 1, i, 1);
        [g.MVARs{j}, g.meanMVAR{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mvar', 1, i, 1);
        % pvarEC(i) no exogenous 
        [g.PVARECs{j}, g.meanPVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'pvarec', 1, j, 1);
        % mpcvarEC(i) auto exogenous 
        [g.MPCVARECs{j}, g.meanMPCVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mpcvarec', 1, i, 1);
        % ppcvarEC(i) no exogenous 
        [g.PPCVARECs{j}, g.meanPPCVAREC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'ppcvarec', 1, j, 1);
        % mpcvarGC(i) no exogenous 
        [g.MPCVARGCs{j}, g.meanMPCVARGC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'mpcvargc', 1, j, 1);
        % ppcvarGC(i) no exogenous 
        [g.PPCVARGCs{j}, g.meanPPCVARGC{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'ppcvargc', 1, j, 1);
        % DLCM(i)-GC linear auto exogenous
        [g.DL2s{j}, g.meanDL2{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlcm', 0, i, 1, []);
        % DLCM(i)-EC linear auto exogenous
        [g.DLW2s{j}, g.meanDLW2{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlw', 0, i, 1, []);
        % DLCM(i)-GC auto exogenous
        [g.DLs{j}, g.meanDL{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlcm', 0, i, 1);
        % DLCM(i)-EC auto exogenous
        [g.DLWs{j}, g.meanDLW{j}, ~] = calculateConnectivity(signals, roiNames, groupName, 'dlw', 0, i, 1);
    end
    % FC no exogenous (pairwise, then exogenous does not have meaning)
    [g.FCs, g.meanFC, ~] = calculateConnectivity(signals, roiNames, groupName, 'fc', 1, j, 0);
end

function statisticalGroupIdentificationByAlgorithms(g, p, roiNames, name1, name2, maxLag)
    global resultsPath;
    global resultsPrefix;

    % plot correlation and cos similarity
    nanx = eye(size(g.meanGC{1},1),size(g.meanGC{1},2));
    nanx(nanx==1) = NaN;
    cosSim = zeros(maxLag*6*2,1);
    figure; bar(cosSim);
    title(['cos similarity between ' name1 ' and ' name2 ' by each algorithm']);

    % normality test
%    g.DLWsNt = calculateAlzNormalityTest(g.DLWs{j}, roiNames, name1, 'dlw');
%    p.DLWsNt = calculateAlzNormalityTest(p.DLWs{j}, roiNames, name2, 'dlw');

    % compalizon test (Wilcoxon, Mann?Whitney U test)
    for j=1:maxLag*2
        [~, GCsUtP{j}, ~] = calculateAlzWilcoxonTest(g.GCs{j}, p.GCs{j}, roiNames, name1, name2, ['gc' num2str(j)]);
        [~, MvarECsUtP{j}, ~] = calculateAlzWilcoxonTest(g.MVARECs{j}, p.MVARECs{j}, roiNames, name1, name2, ['mvarec' num2str(j)]);
        [~, MvarsUtP{j}, ~] = calculateAlzWilcoxonTest(g.MVARs{j}, p.MVARs{j}, roiNames, name1, name2, ['mvar' num2str(j)]);
        [~, MpcvarECsUtP{j}, ~] = calculateAlzWilcoxonTest(g.MPCVARECs{j}, p.MPCVARECs{j}, roiNames, name1, name2, ['mpcvarec' num2str(j)]);
        [~, DL2sUtP{j}, ~] = calculateAlzWilcoxonTest(g.DL2s{j}, p.DL2s{j}, roiNames, name1, name2, ['dlcm_lin' num2str(j)]);
        [~, DLW2sUtP{j}, ~] = calculateAlzWilcoxonTest(g.DLW2s{j}, p.DLW2s{j}, roiNames, name1, name2, ['dlw_lin' num2str(j)]);
        [~, DLsUtP{j}, ~] = calculateAlzWilcoxonTest(g.DLs{j}, p.DLs{j}, roiNames, name1, name2, ['dlcm' num2str(j)]);
        [~, DLWsUtP{j}, ~] = calculateAlzWilcoxonTest(g.DLWs{j}, p.DLWs{j}, roiNames, name1, name2, ['dlw' num2str(j)]);
    end
    [~, FCsUtP, ~] = calculateAlzWilcoxonTest(g.FCs, p.FCs, roiNames, name1, name2, 'fc');

    % using minimum 100 p-value relations. perform 5-fold cross validation.
    topNum = 100;
    sigTh = 2;
    N = 5;

    fcAUC = zeros(1,N);
    gcAUC = zeros(maxLag*2,N);
    mvarecAUC = zeros(maxLag*2,N);
    mvarAUC = zeros(maxLag*2,N);
    mpcvarecAUC = zeros(maxLag*2,N);
    dlAUC = zeros(maxLag*2,N);
    dlwAUC = zeros(maxLag*2,N);
    dl2AUC = zeros(maxLag*2,N);
    dlw2AUC = zeros(maxLag*2,N);
    fcROC = cell(N,2);
    fcACC = cell(N,1);
    for lags=1:maxLag*2
        gcROC{lags} = cell(N,2);
        mvarecROC{lags} = cell(N,2);
        mvarROC{lags} = cell(N,2);
        mpcvarecROC{lags} = cell(N,2);
        dlROC{lags} = cell(N,2);
        dlwROC{lags} = cell(N,2);
        dl2ROC{lags} = cell(N,2);
        dlw2ROC{lags} = cell(N,2);
        gcACC{lags} = cell(N,1);
        mvarecACC{lags} = cell(N,1);
        mvarACC{lags} = cell(N,1);
        mpcvarecACC{lags} = cell(N,1);
        dlACC{lags} = cell(N,1);
        dlwACC{lags} = cell(N,1);
        dl2ACC{lags} = cell(N,1);
        dlw2ACC{lags} = cell(N,1);
    end

    sigCntCN = cell(N,maxLag*2*6);
    sigCntAD = cell(N,maxLag*2*6);
    for k=1:N
        for j=1:maxLag*2
            i = 1;
            % check sigma of healthy subject
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.GCs{j}, p.GCs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, GCsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [gcROC{j}{k,1}, gcROC{j}{k,2}, gcAUC(j,k), gcACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.MVARECs{j}, p.MVARECs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, MvarECsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mvarecROC{j}{k,1}, mvarecROC{j}{k,2}, mvarecAUC(j,k), mvarecACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.MVARs{j}, p.MVARs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, MvarsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mvarROC{j}{k,1}, mvarROC{j}{k,2}, mvarAUC(j,k), mvarACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.MPCVARECs{j}, p.MPCVARECs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, MpcvarECsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mpcvarecROC{j}{k,1}, mpcvarecROC{j}{k,2}, mpcvarecAUC(j,k), mpcvarecACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.DL2s{j}, p.DL2s{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, DL2sUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dl2ROC{j}{k,1}, dl2ROC{j}{k,2}, dl2AUC(j,k), dl2ACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.DLW2s{j}, p.DLW2s{j}, k, N);         % replece g.*s, p.*s
            [B, I, X] = sortAndPairPValues(control, target, DLW2sUtP{j}, topNum);                                  % replace g.p.*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlw2ROC{j}{k,1}, dlw2ROC{j}{k,2}, dlw2AUC(j,k), dlw2ACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.DLs{j}, p.DLs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, DLsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlROC{j}{k,1}, dlROC{j}{k,2}, dlAUC(j,k), dlACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.DLWs{j}, p.DLWs{j}, k, N);         % replece g.*s, p.*s
            [B, I, X] = sortAndPairPValues(control, target, DLWsUtP{j}, topNum);                                  % replace g.p.*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlwROC{j}{k,1}, dlwROC{j}{k,2}, dlwAUC(j,k), dlwACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
        end

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(g.FCs, p.FCs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, FCsUtP, topNum);
        sigCntG1{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntG2{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [fcROC{k,1}, fcROC{k,2}, fcAUC(1,k), fcACC{k}] = calcAlzROCcurve(sigCntG1{k,i}, sigCntG2{k,i}, topNum);
    end

    % save result
    fname = [resultsPath '/' resultsPrefix '-' name1 '-' name2 '-roi' num2str(132) '-result.mat'];
    save(fname, 'cosSim', 'fcAUC','gcAUC','mvarecAUC','mvarAUC','mpcvarecAUC','dlAUC','dlwAUC','dl2AUC','dlw2AUC', ...
        'fcROC','gcROC','mvarecROC','mvarROC','mpcvarecROC','dlROC','dlwROC','dl2ROC','dlw2ROC', ...
        'fcACC','gcACC','mvarecACC','mvarACC','mpcvarecACC','dlACC','dlwACC','dl2ACC','dlw2ACC', ...
        'sigCntCN', 'sigCntAD');

    % show box plot
    AUCs = nan(N,80);
    r = [1:10];
    AUCs(:,r) = gcAUC.'; r=r+10;
    AUCs(:,r) = mvarecAUC.'; r=r+10;
    AUCs(:,r) = mvarAUC.'; r=r+10;
    AUCs(:,r) = mpcvarecAUC.'; r=r+10;
    AUCs(:,r) = dl2AUC.'; r=r+10;
    AUCs(:,r) = dlw2AUC.'; r=r+10;
    AUCs(:,r) = dlAUC.'; r=r+10;
    AUCs(:,r) = dlwAUC.'; r=r+10;
    AUCs(:,r(1)) = fcAUC.';
    figure; boxplot(AUCs);
    title(['AUC box plot : ' name1 ' vs ' name2]);

    % show average ROC curves
    figure; 
    hold on;
    for lags=1:maxLag
        plotAverageROCcurve(gcROC{lags}, N, '-', [0.2,0.5,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mvarecROC{lags}, N, '-', [0.5,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mvarROC{lags}, N, '--', [0.5,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mpcvarecROC{lags}, N, '-.', [0.5,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(dlROC{lags}, N, '-', [0.2,0.2,0.4]+(lags*0.1),1.0);
        plotAverageROCcurve(dl2ROC{lags}, N, '--', [0.2,0.2,0.4]+(lags*0.1),0.4); % linear
        plotAverageROCcurve(dlwROC{lags}, N, '-', [0.2,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(dlw2ROC{lags}, N, '--', [0.2,0.2,0.2]+(lags*0.1),0.4); % linear
    end
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve (without exogenous)']);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')

    % show average ROC curves
    figure; 
    hold on;
    for lags=1:maxLag
        k = lags+5;
        plotAverageROCcurve(gcROC{k}, N, '-', [0.2,0.5,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mvarecROC{k}, N, '-', [0.5,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mvarROC{k}, N, '--', [0.5,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mpcvarecROC{k}, N, '-.', [0.5,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(dlROC{k}, N, '-', [0.2,0.2,0.4]+(lags*0.1),1.0);
        plotAverageROCcurve(dl2ROC{k}, N, '--', [0.2,0.2,0.4]+(lags*0.1),0.4); % linear
        plotAverageROCcurve(dlwROC{k}, N, '-', [0.2,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(dlw2ROC{k}, N, '--', [0.2,0.2,0.2]+(lags*0.1),0.4); % linear
    end
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve (with exogenous)']);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

function [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(orgControl, orgTarget, k, N)
    un = floor(size(orgControl,3) / N);
    st = (k-1)*un+1;
    ed = k*un;
    if k==N, ed = size(orgControl,3); end
    control = orgControl(:,:,st:ed);
    un = floor(size(orgTarget,3) / N);
    st = (k-1)*un+1;
    ed = k*un;
    if k==N, ed = size(orgTarget,3); end
    target = orgTarget(:,:,st:ed);
    if N > 1
        orgTarget(:,:,st:ed) = [];
    end
    meanTarget = nanmean(orgTarget, 3);
    stdTarget = nanstd(orgTarget, 1, 3);
    meanControl = nanmean(orgControl, 3);
end

function [x, y, auc] = invertROCcurve(inx, iny)
    y = inx;
    x = iny;
    auc = trapz(x, y);
end

function [x, y, auc, accuracy] = calcAlzROCcurve(control, target, start)
    x = [0]; y = [0]; % start from (0,0)
    accuracy = nan(start+1,1);
    tpmax = length(control);
    fpmax = length(target);
    for i=start:-1:0
        tp = length(find(control>=i));
        fp = length(find(target>=i));
        tn = fpmax - fp;
        x = [x fp/fpmax];
        y = [y tp/tpmax];
        accuracy(i+1) = (tp + tn) / (tpmax + fpmax);
    end
    auc = trapz(x, y);
end

function [B, I, X] = sortAndPairPValues(control, target, utestP2, topNum)
    ROINUM = size(control,1);
    [B, I] = sort(utestP2(:));
    X = [];
    for k=1:topNum
        i = floor(mod(I(k)-1,ROINUM) + 1);
        j = floor((I(k)-1)/ROINUM + 1);
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
        X = [X, x2, y2];
    end
end

function sigCount = calcAlzSigmaSubjects(weights, meanWeights, stdWeights, meanControl, I, topNum, sigTh)
    ROINUM = size(weights,1);
    subjectNum = size(weights,3);
    X = [];
    sigCount = [];
    isControlBig = meanControl - meanWeights;
    for n=1:subjectNum
        w2 = weights(:,:,n);
        w2 = w2 - meanWeights;
        sig = abs(w2 ./ stdWeights);
%        sig = w2 ./ stdWeights;
        s = nan(topNum, 1);
        for k=1:topNum
            i = floor(mod(I(k)-1,ROINUM) + 1);
            j = floor((I(k)-1)/ROINUM + 1);
            s(k) = sig(i, j);
%{
            s2 = sig(i, j);
            if s2 < 0 && isControlBig(i, j) > 0
                s2 = 0;
            elseif s2 > 0 && isControlBig(i, j) < 0
                s2 = 0;
            end
            s(k) = abs(s2);
%}
        end
        X = [X, s];
        sigCount = [sigCount, length(find(s>=sigTh))];
    end
%    figure; boxplot(X);
%    figure; bar(sigCount);
end

function [normalities, normalitiesP] = calculateAlzNormalityTest(weights, roiNames, group, algorithm)
    % constant value
    ROINUM = size(weights,1);

    global resultsPath;
    global resultsPrefix;
    outfName = [resultsPath '/' resultsPrefix '-' algorithm '-' group '-roi' num2str(ROINUM) '-normality.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        normalities = nan(ROINUM, ROINUM);
        normalitiesP = nan(ROINUM, ROINUM);
        for i=1:ROINUM
            for j=1:ROINUM
                if i==j, continue; end
                x = squeeze(weights(i,j,:));
                [h, p] = lillietest(x);
                normalities(i,j) = 1 - h;
                normalitiesP(i,j) = p;
            end
        end
        save(outfName, 'normalities', 'normalitiesP', 'roiNames');
    end

    load('test/colormap.mat')
    % show normality test result
    figure; 
    colormap(hvalmap);
    clims = [0,1];
    imagesc(normalities,clims);
    daspect([1 1 1]);
    title([group '-' algorithm ' normality test result']);
    colorbar;
    % normality test p values
    normalitiesP(isnan(normalitiesP)) = 1;
    figure;
    colormap(pvalmap);
    clims = [0,0.5];
    imagesc(normalitiesP,clims);
    daspect([1 1 1]);
    title([group '-' algorithm ' normality test p values']);
    colorbar;
end


function [utestH, utestP, utestP2] = calculateAlzWilcoxonTest(control, target, roiNames, controlGroup, targetGroup, algorithm)
    % constant value
    ROINUM = size(control,1);

    global resultsPath;
    global resultsPrefix;
    outfName = [resultsPath '/' resultsPrefix '-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROINUM) '-utest.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        utestH = nan(ROINUM, ROINUM);
        utestP = nan(ROINUM, ROINUM);
        utestP2 = nan(ROINUM, ROINUM);
        for i=1:ROINUM
            for j=1:ROINUM
                if i==j, continue; end
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
                [p, h] = ranksum(x2,y2);
%                [p, h] = signrank(x2,y2);
                utestH(i,j) = h;
                utestP(i,j) = p;
                if h > 0 && nanmean(x) > nanmean(y)
                    utestP2(i,j) = p;
                end
            end
        end
        save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames');
    end
    % counting by source region and target region
    countSource = nansum(utestH,1);
    countTarget = nansum(utestH,2);
    save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames', 'countSource', 'countTarget');

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
