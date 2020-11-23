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

    % calculate zi & zi\{j}
    [cnDLWs, cnDLWnss, meanCnDLWns, stdCnDLWns, cnInSignals, cnInControls] = calculateDistributions(cnSignals, roiNames, 'cn', 'dlw');
    [adDLWs, adDLWnss, meanAdDLWns, stdAdDLWns, ~, ~] = calculateDistributions(adSignals, roiNames, 'ad', 'dlw');

    % transform healthy zi,zi\{j} to ad's (type 2)
    meanCnDLWns3 = repmat(meanCnDLWns,[1 1 size(cnDLWnss,3)]);
    stdCnDLWns3 = repmat(stdCnDLWns,[1 1 size(cnDLWnss,3)]);
    cnDLWsig = (cnDLWnss - meanCnDLWns3) ./ stdCnDLWns3;
    meanAdDLWns3 = repmat(meanAdDLWns,[1 1 size(cnDLWnss,3)]);
    stdAdDLWns3 = repmat(stdAdDLWns,[1 1 size(cnDLWnss,3)]);
    vadDLWnss = meanAdDLWns3 + cnDLWsig .* stdAdDLWns3;

    % calculate virtual AD ECcnDLWnss (type 2)
    vadDLWs = vadDLWnss(:,2:end,:);
    vadZi = repmat(vadDLWnss(:,1,:),[1 size(vadDLWs,2) 1]);
    vadDLWs = abs(vadZi - vadDLWs);

    % retraining DLCM network (type 2)
    [vad2DLWs, meanVad2DLWns, stdVad2DLWns] = retrainDLCMAndEC(vadDLWnss, roiNames, 'vadns');

    % generate virtual ad signals (type 3)
    vadSignals = calculateVirtualADSignals3(cnSignals, roiNames, cnDLWs, adDLWs, 'dlw');
    [vad3DLs, ~, ~] = calculateConnectivity(vadSignals, roiNames, 'vad3', 'dlcm');
    [vad3DLWs, ~, ~] = calculateConnectivity(vadSignals, roiNames, 'vad3', 'dlw');
    [~, vad3DLWnss, meanVad3DLWns, stdVad3DLWns, ~, ~] = calculateDistributions(vadSignals, roiNames, 'vad3', 'dlw');

    % generate virtual ad signals (type 4)
    [vad4Signals, vad4DLWs, vad4DLWnss] = calculateVirtualADSignals4(cnSignals, adSignals, roiNames, cnInSignals, cnInControls, 'vad4');
    [vad5Signals, vad5DLWs, vad5DLWnss] = calculateVirtualADSignals4(vad4Signals, adSignals, roiNames, cnInSignals, cnInControls, 'vad5');
    [vad6Signals, vad6DLWs, vad6DLWnss] = calculateVirtualADSignals4(vad5Signals, adSignals, roiNames, cnInSignals, cnInControls, 'vad6');

    % plot correlation and cos similarity
    algNum = 8;
    meanCnDLW = nanmean(cnDLWs,3);
    meanAdDLW = nanmean(adDLWs,3);
    meanVadDLW = nanmean(vadDLWs,3);
    meanVad2DLW = nanmean(vad2DLWs,3);
    meanVad3DLW = nanmean(vad3DLWs,3);
    meanVad4DLW = nanmean(vad4DLWs,3);
%    nanx = eye(size(meanCnDLW,1),size(meanCnDLW,2));
%    nanx(nanx==1) = NaN;
%    figure; cnadDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanAdDLW);
%    figure; cnvadDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanVadDLW);
%    figure; advadDLWr = plotTwoSignalsCorrelation(meanAdDLW, meanVadDLW);
%    figure; advadDLWr2 = plotTwoSignalsCorrelation(meanAdDLW, meanVad2DLW);
%    figure; advadDLWr3 = plotTwoSignalsCorrelation(meanAdDLW, meanVad3DLW);
    figure; advadDLWr4 = plotTwoSignalsCorrelation(meanAdDLW, meanVad4DLW);
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanAdDLW);
    cosSim(2) = getCosSimilarity(meanCnDLW, meanVadDLW);
    cosSim(3) = getCosSimilarity(meanAdDLW, meanVadDLW);
    cosSim(4) = getCosSimilarity(meanAdDLW, meanVad2DLW);
    cosSim(5) = getCosSimilarity(meanCnDLW, meanVad3DLW);
    cosSim(6) = getCosSimilarity(meanAdDLW, meanVad3DLW);
    cosSim(7) = getCosSimilarity(meanCnDLW, meanVad4DLW);
    cosSim(8) = getCosSimilarity(meanAdDLW, meanVad4DLW);
    X = categorical({'cn-ad','cn-vad','ad-vad','ad-vad2','cn-vad3','ad-vad3','cn-vad4','ad-vad4'});
    figure; bar(X, cosSim);
    title('cos similarity between CN and AD by each algorithm');

    % normality test
%{
    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cnec', 'dlw');
    adDLWsNt = calculateAlzNormalityTest(adDLWs, roiNames, 'adec', 'dlw');
    vadDLWsNt = calculateAlzNormalityTest(vadDLWs, roiNames, 'vadec', 'dlw');
    cnnsDLWsNt = calculateAlzNormalityTest(cnDLWnss, roiNames, 'cnns', 'dlw');
    adnsDLWsNt = calculateAlzNormalityTest(adDLWnss, roiNames, 'adns', 'dlw');
    vadnsDLWsNt = calculateAlzNormalityTest(vadDLWnss, roiNames, 'vadns', 'dlw');
    vad3nsDLWsNt = calculateAlzNormalityTest(vad3DLWnss, roiNames, 'vad3ns', 'dlw');
    vad4nsDLWsNt = calculateAlzNormalityTest(vad4DLWnss, roiNames, 'vad4ns', 'dlw');
%}
    % compalizon test (Wilcoxon, Mann?Whitney U test)
    [cnadDLWsUt, cnadDLWsUtP, cnadDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, adDLWs, roiNames, 'cnec', 'adec', 'dlw');
    [cnvadDLWsUt, cnvadDLWsUtP, cnvadDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vadDLWs, roiNames, 'cnec', 'vadec', 'dlw');
    [advadDLWsUt, advadDLWsUtP, advadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vadDLWs, roiNames, 'adec', 'vadec', 'dlw');
    [advad2DLWsUt, advad2DLWsUtP, advad2DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad2DLWs, roiNames, 'adec', 'vad2ec', 'dlw');
    [advad3DLWsUt, advad3DLWsUtP, advad3DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad3DLWs, roiNames, 'adec', 'vad3ec', 'dlw');
    [cnadDLWnssUt, cnadDLWnssUtP, cnadDLWnssUtP2] = calculateAlzWilcoxonTest(cnDLWnss, adDLWnss, roiNames, 'cnns', 'adns', 'dlw');
    [cnvadDLWnssUt, cnvadDLWnssUtP, cnvadDLWnssUtP2] = calculateAlzWilcoxonTest(cnDLWnss, vadDLWnss, roiNames, 'cnns', 'vadns', 'dlw');
    [advadDLWnssUt, advadDLWnssUtP, advadDLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vadDLWnss, roiNames, 'adns', 'vadns', 'dlw');
    [advad3DLWnssUt, advad3DLWnssUtP, advad3DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad3DLWnss, roiNames, 'adns', 'vad3ns', 'dlw');
    [advad4DLWnssUt, advad4DLWnssUtP, advad4DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad4DLWnss, roiNames, 'adns', 'vad4ns', 'dlw');
    
%{
    % using minimum 100 p-value relations. perform 5-fold cross validation.
    topNum = 100;
    sigTh = 2;
    N = 1;

    dlwAUC = zeros(1,N);
    dlwvAUC = zeros(1,N);
    dlwROC = cell(N,2);
    dlwvROC = cell(N,2);

    sigCntCN = cell(N,algNum);
    sigCntAD = cell(N,algNum);
    for k=1:N
        i = 1;
        % check sigma of healthy subject
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs, adDLWs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadDLWsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

        % check sigma of healthy subject
        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs, vadDLWs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnvadDLWsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlwvROC{k,1}, dlwvROC{k,2}, dlwvAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
    end
    figure; boxplot(X);

    % save result
    fname = ['results/adsim-cn-ad-roi' num2str(132) '-result.mat'];
    save(fname, 'cosSim', 'dlwAUC','dlwvAUC','dlwROC','dlwvROC', 'sigCntCN', 'sigCntAD');
    mean(dlwAUC) % show result AUC
    mean(dlwvAUC) % show result AUC

    % plot ROC curve
    figure;
    hold on;
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % TODO:
    plotErrorROCcurve(dlwvROC, N, [0.6,0.2,0.2]); % TODO:
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.0);
    plotAverageROCcurve(dlwvROC, N, '-', [0.6,0.2,0.2],1.0);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve']);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
%}
end

function [weights, meanWeights, stdWeights] = retrainDLCMAndEC(signals, roiNames, group)
    ROWNUM = size(signals,1);
    COLNUM = size(signals,2);
    sbjNum = size(signals,3);
    weights = zeros(ROWNUM, ROWNUM, sbjNum);

    outfName = ['results/adsim-retrain-' group '-roi' num2str(ROWNUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        for i=1:sbjNum
            dlcmName = ['results/adsim-dlcm-' group '-roi' num2str(ROWNUM) '-net' num2str(i) '.mat'];
            if exist(dlcmName, 'file')
                load(dlcmName);
            else
                nodeInputOrg = ones(ROWNUM*2, COLNUM);
                nodeInputOrg(1:ROWNUM,2:end) = nodeInputOrg(1:ROWNUM, 2:end) - eye(ROWNUM);

                origName = ['results/ad-dlcm-cn-roi' num2str(ROWNUM) '-net' num2str(i) '.mat'];
                load(origName);
                if exist('m','var')
                    c = m;
                end

                % training DLCM network
                maxEpochs = 30;
                options = trainingOptions('adam', ...
                    'ExecutionEnvironment','cpu', ...
                    'MaxEpochs',maxEpochs, ...
                    'MiniBatchSize',ROWNUM, ...
                    'Shuffle','every-epoch', ...
                    'GradientThreshold',5,...
                    'L2Regularization',0.05, ...
                    'Verbose',false);
    %                'Plots','training-progress');

                disp('start training');
                for j=1:ROWNUM
                    disp(['virtual alzheimer training node ' num2str(i) '-' num2str(j)]);
                    nodeTeach = signals(j,1:end,i);
                    nodeInput = nodeInputOrg;
                    if ~isempty(inControl)
                        filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
                        nodeInput(ROWNUM+1:end,:) = nodeInput(ROWNUM+1:end,:) .* filter;
                    end
                    nodeTeach(:,j+1) = [];
                    nodeInput(:,j+1) = [];
                    [netDLCM.nodeNetwork{j}, netDLCM.trainInfo{j}] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{j}, options);
                end

                save(dlcmName, 'netDLCM', 'si', 'inSignal', 'inControl', 'mat', 'sig', 'c', 'maxsi', 'minsi');
            end

            % recalculate EC
            weights(:,:,i) = calcDlcmEC(netDLCM, [], inControl);
        end
        save(outfName, 'weights', 'roiNames');
    end
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);
end

function vadSignals = calculateVirtualADSignals3(signals, roiNames, cnECs, adECs, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    vadSignals = {};
    
    outfName = ['results/adsim-signal-vad-roi' num2str(ROINUM) '-' algorithm '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end
    
    cnOutECs = nansum(nanmean(cnECs, 3), 1);
    adOutECs = nansum(nanmean(adECs, 3), 1);
    transRate = adOutECs ./ cnOutECs;

    for i=1:length(signals)
        sbjSignals = signals{i};
        for j=1:ROINUM
            sbjSignals(j,:) = sbjSignals(j,:) * transRate(1,j);
        end
        vadSignals{end+1} = sbjSignals;
    end
    save(outfName, 'vadSignals', 'roiNames');
end

function [vadSignals, vadDLWs, vadDLWnss] = calculateVirtualADSignals4(cnSignals, adSignals, roiNames, cnInSignals, cnInControls, group)
    ROWNUM = size(cnSignals{1},1);
    sigLen = size(cnSignals{1},2);
    cnNum = length(cnSignals);
    adNum = length(adSignals);
    outfName = ['results/adsim-all-vad-roi' num2str(ROWNUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        allCnSignals = nan(ROWNUM, sigLen, cnNum);
        allVadSignals = nan(ROWNUM, sigLen, cnNum, adNum);
        sig = nan(cnNum);
        c = nan(cnNum);
        maxsi = nan(cnNum);
        minsi = nan(cnNum);
        for k=1:adNum
            disp(['generate virtual ad signals (type 4) : ' num2str(k)]);
            origName = ['results/ad-dlcm-ad-roi' num2str(ROWNUM) '-net' num2str(k) '.mat'];
            load(origName);
            for i=1:cnNum
                [allCnSignals(:,:,i), sig(i), c(i), maxsi(i), minsi(i)] = convert2SigmoidSignal(cnSignals{i});
                [S, time, mae, maeerr] = predictDlcmNetwork(allCnSignals(:,:,i), cnInSignals(:,:,i), [], cnInControls(:,:,1), netDLCM);
                allVadSignals(:,:,i,k) = S;
            end
        end
        save(outfName, 'allCnSignals', 'cnInSignals', 'sig', 'c', 'maxsi', 'minsi', 'allVadSignals');
    end
    % get mean of AD DLCM generated signals (type 4)
    meanVadSignals = nanmean(allVadSignals, 4);
    vadSignals = {};
    for i=1:cnNum
        vadSignals{end+1} = convert2InvSigmoidSignal(meanVadSignals(:,:,i), sig(i), c(i), maxsi(i), minsi(i));
%        plot(vad4Signals{end});
    end
    [vad4DLs, ~, ~] = calculateConnectivity(vadSignals, roiNames, group, 'dlcm');
    [vadDLWs, ~, ~] = calculateConnectivity(vadSignals, roiNames, group, 'dlw');
    [~, vadDLWnss, meanVad4DLWns, stdVad4DLWns, ~, ~] = calculateDistributions(vadSignals, roiNames, group, 'dlw');
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

function [x, y, auc] = calcAlzROCcurve(control, target, start)
    x = [0]; y = [0]; % start from (0,0)
    tpmax = length(control);
    fpmax = length(target);
    for i=start:-1:0
        tp = length(find(control>=i));
        fp = length(find(target>=i));
        x = [x fp/fpmax];
        y = [y tp/tpmax];
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

function [ECs, nodeSignals, meanSignals, stdSignals, inSignals, inControls] = calculateDistributions(signals, roiNames, group, algorithm)
    % constant value
    ROINUM = size(signals{1},1);

    outfName = ['results/adsim-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        ECs = zeros(ROINUM, ROINUM, length(signals));
        nodeSignals = zeros(ROINUM, ROINUM+1, length(signals));
        inSignals = zeros(ROINUM, size(signals{1},2), length(signals));
        inControls = zeros(ROINUM, ROINUM, length(signals));
        for i=1:length(signals)
            switch(algorithm)
            case 'dlw'
                dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                load(dlcmName);
                [ec, ecNS] = calcDlcmEC(netDLCM, [], inControl);
            end
            ECs(:,:,i) = ec;
            nodeSignals(:,:,i) = ecNS;
            inSignals(:,:,i) = inSignal;
            inControls(:,:,i) = inControl;
        end
        save(outfName, 'ECs', 'nodeSignals', 'roiNames', 'inSignals', 'inControls');
    end
    meanSignals = nanmean(nodeSignals, 3);
    stdSignals = nanstd(nodeSignals, 1, 3);
    save(outfName, 'ECs', 'nodeSignals', 'meanSignals', 'stdSignals', 'roiNames', 'inSignals', 'inControls');
end

function [normalities, normalitiesP] = calculateAlzNormalityTest(ECs, roiNames, group, algorithm)
    % constant value
    ROWNUM = size(ECs,1);
    COLNUM = size(ECs,2);

    outfName = ['results/adsim-' algorithm '-' group '-roi' num2str(ROWNUM) '-normality.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        normalities = nan(ROWNUM, COLNUM);
        normalitiesP = nan(ROWNUM, COLNUM);
        for i=1:ROWNUM
            for j=1:COLNUM
                if isnan(ECs(i,j,1)), continue; end
                x = squeeze(ECs(i,j,:));
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
    ROWNUM = size(control,1);
    COLNUM = size(control,2);

    outfName = ['results/adsim-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROWNUM) '-utest.mat'];
    if exist(outfName, 'file')
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
                [p, h] = ranksum(x2,y2);
                %[p, h] = signrank(x2,y2);
                %[h, p] = ttest2(x2,y2);
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
