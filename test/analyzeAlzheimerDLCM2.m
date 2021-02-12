function analyzeAlzheimerDLCM2
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};
    pathesMCI = {'ADNI2_65-75_F_MCI_nii', 'ADNI2_65-75_M_MCI_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');
    [mciSignals] = connData2signalsFile(base, pathesMCI, 'mci');

    maxLag = 5;
    % calculate connectivity
    for j=1:maxLag
        % mvGC(i) no exogenous 
        [cnGCs{j}, meanCNGC{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'gc', 1, j, 0);
        [adGCs{j}, meanADGC{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'gc', 1, j, 0);
        % mvarEC(i) no exogenous 
        [cnMVARECs{j}, meanCNMVAREC{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mvarec', 1, j, 0);
        [adMVARECs{j}, meanADMVAREC{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'mvarec', 1, j, 0);
        % DLCM(i)-GC linear no exogenous
        [cnDL2s{j}, meanCNDL2{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 0, j, 0, []);
        [adDL2s{j}, meanADDL2{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 0, j, 0, []);
        % DLCM(i)-EC linear no exogenous
        [cnDLW2s{j}, meanCNDLW2{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlw', 0, j, 0, []);
        [adDLW2s{j}, meanADDLW2{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlw', 0, j, 0, []);
        % DLCM(i)-GC no exogenous
        [cnDLs{j}, meanCNDL{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 0, j, 0);
        [adDLs{j}, meanADDL{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 0, j, 0);
        % DLCM(i)-EC no exogenous
        [cnDLWs{j}, meanCNDLW{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlw', 0, j, 0);
        [adDLWs{j}, meanADDLW{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlw', 0, j, 0);
    end

    for i=1:maxLag
        j = i+maxLag;
        % mvGC(i) auto exogenous 
        [cnGCs{j}, meanCNGC{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'gc', 1, i, 1);
        [adGCs{j}, meanADGC{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'gc', 1, i, 1);
        % mvarEC(i) auto exogenous 
        [cnMVARECs{j}, meanCNMVAREC{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mvarec', 1, i, 1);
        [adMVARECs{j}, meanADMVAREC{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'mvarec', 1, i, 1);
        % DLCM(i)-GC linear auto exogenous
        [cnDL2s{j}, meanCNDL2{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 0, i, 1, []);
        [adDL2s{j}, meanADDL2{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 0, i, 1, []);
        % DLCM(i)-EC linear auto exogenous
        [cnDLW2s{j}, meanCNDLW2{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlw', 0, i, 1, []);
        [adDLW2s{j}, meanADDLW2{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlw', 0, i, 1, []);
        % DLCM(i)-GC auto exogenous
        [cnDLs{j}, meanCNDL{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 0, i, 1);
        [adDLs{j}, meanADDL{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 0, i, 1);
        % DLCM(i)-EC auto exogenous
        [cnDLWs{j}, meanCNDLW{j}, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlw', 0, i, 1);
        [adDLWs{j}, meanADDLW{j}, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlw', 0, i, 1);
    end
    
    % plot correlation and cos similarity
    nanx = eye(size(meanCNGC{1},1),size(meanCNGC{1},2));
    nanx(nanx==1) = NaN;
    cosSim = zeros(maxLag*6*2,1);
    figure; bar(cosSim);
    title('cos similarity between CN and AD by each algorithm');

    % normality test
%    cnDLWsNt = calculateAlzNormalityTest(cnDLWs{j}, roiNames, 'cn', 'dlw');
%    adDLWsNt = calculateAlzNormalityTest(adDLWs{j}, roiNames, 'ad', 'dlw');
%    mciDLWsNt = calculateAlzNormalityTest(mciDLWs{j}, roiNames, 'mci', 'dlw');

    % compalizon test (Wilcoxon, Mann?Whitney U test)
    for j=1:maxLag*2
        [~, cnadGCsUtP{j}, ~] = calculateAlzWilcoxonTest(cnGCs{j}, adGCs{j}, roiNames, 'cn', 'ad', ['gc' num2str(j)]);
        [~, cnadMvarECsUtP{j}, ~] = calculateAlzWilcoxonTest(cnMVARECs{j}, adMVARECs{j}, roiNames, 'cn', 'ad', ['mvarec' num2str(j)]);
        [~, cnadDL2sUtP{j}, ~] = calculateAlzWilcoxonTest(cnDL2s{j}, adDL2s{j}, roiNames, 'cn', 'ad', ['dlcm_lin' num2str(j)]);
        [~, cnadDLW2sUtP{j}, ~] = calculateAlzWilcoxonTest(cnDLW2s{j}, adDLW2s{j}, roiNames, 'cn', 'ad', ['dlw_lin' num2str(j)]);
        [~, cnadDLsUtP{j}, ~] = calculateAlzWilcoxonTest(cnDLs{j}, adDLs{j}, roiNames, 'cn', 'ad', ['dlcm' num2str(j)]);
        [~, cnadDLWsUtP{j}, ~] = calculateAlzWilcoxonTest(cnDLWs{j}, adDLWs{j}, roiNames, 'cn', 'ad', ['dlw' num2str(j)]);
    end

    % using minimum 100 p-value relations. perform 5-fold cross validation.
    topNum = 100;
    sigTh = 2;
    N = 5;

    gcAUC = zeros(maxLag*2,N);
    mvarecAUC = zeros(maxLag*2,N);
    dlAUC = zeros(maxLag*2,N);
    dlwAUC = zeros(maxLag*2,N);
    dl2AUC = zeros(maxLag*2,N);
    dlw2AUC = zeros(maxLag*2,N);
    for lags=1:maxLag*2
        gcROC{lags} = cell(N,2);
        mvarecROC{lags} = cell(N,2);
        dlROC{lags} = cell(N,2);
        dlwROC{lags} = cell(N,2);
        dl2ROC{lags} = cell(N,2);
        dlw2ROC{lags} = cell(N,2);
        gcACC{lags} = cell(N,1);
        mvarecACC{lags} = cell(N,1);
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
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnGCs{j}, adGCs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, cnadGCsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [gcROC{j}{k,1}, gcROC{j}{k,2}, gcAUC(j,k), gcACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMVARECs{j}, adMVARECs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, cnadMvarECsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mvarecROC{j}{k,1}, mvarecROC{j}{k,2}, mvarecAUC(j,k), mvarecACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDL2s{j}, adDL2s{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, cnadDL2sUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dl2ROC{j}{k,1}, dl2ROC{j}{k,2}, dl2AUC(j,k), dl2ACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLW2s{j}, adDLW2s{j}, k, N);         % replece cn*s, ad*s
            [B, I, X] = sortAndPairPValues(control, target, cnadDLW2sUtP{j}, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlw2ROC{j}{k,1}, dlw2ROC{j}{k,2}, dlw2AUC(j,k), dlw2ACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLs{j}, adDLs{j}, k, N);
            [B, I, X] = sortAndPairPValues(control, target, cnadDLsUtP{j}, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlROC{j}{k,1}, dlROC{j}{k,2}, dlAUC(j,k), dlACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs{j}, adDLWs{j}, k, N);         % replece cn*s, ad*s
            [B, I, X] = sortAndPairPValues(control, target, cnadDLWsUtP{j}, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlwROC{j}{k,1}, dlwROC{j}{k,2}, dlwAUC(j,k), dlwACC{j}{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
        end
    end

    % save result
    fname = ['results/ad-cn-ad-roi' num2str(132) '-result.mat'];
    save(fname, 'cosSim', 'gcAUC','mvarecAUC','dlAUC','dlwAUC','dl2AUC','dlw2AUC', ...
        'gcROC','mvarecROC','dlROC','dlwROC','dl2ROC','dlw2ROC', ...
        'gcACC','mvarecACC','dlACC','dlwACC','dl2ACC','dlw2ACC', ...
        'sigCntCN', 'sigCntAD');

    % show box plot
    AUCs = nan(N,60);
    r = [1:10];
    AUCs(:,r) = gcAUC.'; r=r+10;
    AUCs(:,r) = mvarecAUC.'; r=r+10;
    AUCs(:,r) = dl2AUC.'; r=r+10;
    AUCs(:,r) = dlw2AUC.'; r=r+10;
    AUCs(:,r) = dlAUC.'; r=r+10;
    AUCs(:,r) = dlwAUC.'; r=r+10;
    figure; boxplot(AUCs);
    title('AUC box plot idx');

    % show average ROC curve of DCM
    figure; 
    hold on;
    for lags=1:maxLag
        plotAverageROCcurve(gcROC{lags}, N, '--', [0.2,0.5,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mvarecROC{lags}, N, '-', [0.5,0.2,0.2]+(lags*0.1),1.0);
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

    % show average ROC curve of DCM
    figure; 
    hold on;
    for lags=1:maxLag
        k = lags+5;
        plotAverageROCcurve(gcROC{k}, N, '--', [0.2,0.5,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(mvarecROC{k}, N, '-', [0.5,0.2,0.2]+(lags*0.1),1.0);
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

    outfName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-normality.mat'];
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

    outfName = ['results/ad-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROINUM) '-utest.mat'];
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
