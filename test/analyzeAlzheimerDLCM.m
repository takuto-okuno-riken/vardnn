function analyzeAlzheimerDLCM
    % fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-70_F_CN_nii', 'ADNI2_65-70_M_CN_nii', 'ADNI2_70-75_F_CN_nii', 'ADNI2_70-75_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};
    pathesMCI = {'ADNI2_65-75_F_MCI_nii', 'ADNI2_65-75_M_MCI_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');
    [mciSignals] = connData2signalsFile(base, pathesMCI, 'mci');

    % calculate connectivity
    [cnFCs, meanCNFC, stdCNFC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc');
    [adFCs, meanADFC, stdADFC] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc');
    [mciFCs, meanMCIFC, stdMCIFC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'fc');

    [cnGCs, meanCNGC, stdCNGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'gc');
    [adGCs, meanADGC, stdADGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'gc');
    [mciGCs, meanMCIGC, stdMCIGC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'gc');

    [cnTEs, meanCNTE, stdCNTE] = calculateConnectivity(cnSignals, roiNames, 'cn', 'te');
    [adTEs, meanADTE, stdADTE] = calculateConnectivity(adSignals, roiNames, 'ad', 'te');
    [mciTEs, meanMCITE, stdMCITE] = calculateConnectivity(mciSignals, roiNames, 'mci', 'te');

    [cnDLs, meanCNDL, stdCNDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm');
    [adDLs, meanADDL, stdADDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm');
    [mciDLs, meanMCIDL, stdMCIDL] = calculateConnectivity(mciSignals, roiNames, 'mci', 'dlcm');
    
    % normality test
%{
    cnFCsNt = calculateNormalityTest(cnFCs, roiNames, 'cn', 'fc');
    adFCsNt = calculateNormalityTest(adFCs, roiNames, 'ad', 'fc');
    mciFCsNt = calculateNormalityTest(mciFCs, roiNames, 'mci', 'fc');

    cnGCsNt = calculateNormalityTest(cnGCs, roiNames, 'cn', 'gc');
    adGCsNt = calculateNormalityTest(adGCs, roiNames, 'ad', 'gc');
    mciGCsNt = calculateNormalityTest(mciGCs, roiNames, 'mci', 'gc');

    cnTEsNt = calculateNormalityTest(cnTEs, roiNames, 'cn', 'te');
    adTEsNt = calculateNormalityTest(adTEs, roiNames, 'ad', 'te');
    mciTEsNt = calculateNormalityTest(mciTEs, roiNames, 'mci', 'te');

    cnDLsNt = calculateNormalityTest(cnDLs, roiNames, 'cn', 'dlcm');
    adDLsNt = calculateNormalityTest(adDLs, roiNames, 'ad', 'dlcm');
    mciDLsNt = calculateNormalityTest(mciDLs, roiNames, 'mci', 'dlcm');
%}
    % compalizon test (Wilcoxon, Mann?Whitney U test)
    [cnadFCsUt, cnadFCsUtP, cnadFCsUtP2] = calculateWilcoxonTest(cnFCs, adFCs, roiNames, 'cn', 'ad', 'fc');
    [cnadGCsUt, cnadGCsUtP, cnadGCsUtP2] = calculateWilcoxonTest(cnGCs, adGCs, roiNames, 'cn', 'ad', 'gc');
    [cnadTEsUt, cnadTEsUtP, cnadTEsUtP2] = calculateWilcoxonTest(cnTEs, adTEs, roiNames, 'cn', 'ad', 'te');
    [cnadDLsUt, cnadDLsUtP, cnadDLsUtP2] = calculateWilcoxonTest(cnDLs, adDLs, roiNames, 'cn', 'ad', 'dlcm');
    
    % box plot of minimum p-value relation
    topNum = 50;
    sigTh = 2;

    % check sigma of healthy subject
    [B, I] = boxPlotPValues(cnFCs, adFCs, cnadFCsUtP2, topNum);
    sigCntCNFC = calcSigmaSubjects(cnFCs, meanADFC, stdADFC, I, topNum, sigTh);
    sigCntADFC = calcSigmaSubjects(adFCs, meanADFC, stdADFC, I, topNum, sigTh);
    [B, I] = boxPlotPValues(cnGCs, adGCs, cnadGCsUtP2, topNum);
    sigCntCNGC = calcSigmaSubjects(cnGCs, meanADGC, stdADGC, I, topNum, sigTh);
    sigCntADGC = calcSigmaSubjects(adGCs, meanADGC, stdADGC, I, topNum, sigTh);
    [B, I] = boxPlotPValues(cnTEs, adTEs, cnadTEsUtP2, topNum);
    sigCntCNTE = calcSigmaSubjects(cnTEs, meanADTE, stdADTE, I, topNum, sigTh);
    sigCntADTE = calcSigmaSubjects(adTEs, meanADTE, stdADTE, I, topNum, sigTh);
    [B, I] = boxPlotPValues(cnDLs, adDLs, cnadDLsUtP2, topNum);
    sigCntCNDL = calcSigmaSubjects(cnDLs, meanADDL, stdADDL, I, topNum, sigTh);
    sigCntADDL = calcSigmaSubjects(adDLs, meanADDL, stdADDL, I, topNum, sigTh);
    
    % calculate ROC curve
    figure;
    aucFC = plotAlzROCcurve(sigCntCNFC, sigCntADFC, '-', [0.8,0.2,0.2],0.5)
    aucGC = plotAlzROCcurve(sigCntCNGC, sigCntADGC, '-', [0.1,0.8,0.1],0.5)
    aucTE = plotAlzROCcurve(sigCntCNTE, sigCntADTE, '--', [0.2,0.5,0.7],0.5)
    aucDL = plotAlzROCcurve(sigCntCNDL, sigCntADDL, '-', [0.2,0.2,0.2],1.2)
end

function auc = plotAlzROCcurve(control, target, line, col, width)
    x = [0]; y = [0]; % start from (0,0)
    st = max([control(:); target(:)]);
    tpmax = length(control);
    fpmax = length(target);
    for i=st:-1:0
        tp = length(find(control>=i));
        fp = length(find(target>=i));
        x = [x fp/fpmax];
        y = [y tp/tpmax];
    end
    auc = trapz(x, y);

    % plot ROC curve
    hold on;
    plot(x, y, line,'Color',col,'LineWidth',width);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;    
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title('ROC curve');
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

function [B, I] = boxPlotPValues(control, target, utestP2, topNum)
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
    figure;
    boxplot(X);
end

function sigCount = calcSigmaSubjects(weights, meanWeights, stdWeights, I, topNum, sigTh)
    ROINUM = size(weights,1);
    subjectNum = size(weights,3);
    X = [];
    sigCount = [];
    for n=1:subjectNum
        w2 = weights(:,:,n);
        w2 = w2 - meanWeights;
        sig = abs(w2 ./ stdWeights);
        s = nan(topNum, 1);
        for k=1:topNum
            i = floor(mod(I(k)-1,ROINUM) + 1);
            j = floor((I(k)-1)/ROINUM + 1);
            s(k) = sig(i, j);
        end
        X = [X, s];
        sigCount = [sigCount, length(find(s>=sigTh))];
    end
    figure; boxplot(X);
    figure; bar(sigCount);
end

function [signals, roiNames] = connData2signalsFile(base, pathes, group)
    % constant value
    ROINUM = 132;
    START = 4;
    ReLOCATE = [1	3	5	7	9	11	13	15	17	19	21	23	25	27	29	31	33	35	37	39	41	43	45	47	50	53	58	60	62	64	66	68	70	72	74	76	78	80	82	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	49	52	55	56	57	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	51	54	59	61	63	65	67	69	71	73	75	77	79	81	83	85	87	89	91	93	95	97	99	101	103	105	107	109	111	113	115	117	119	121	123	125	126	127	128	129	130	131	132];


    % init values
    signals = {};
    roiNames2 = cell(1,ROINUM);
    idx = 0;
    
    outfName = ['data/ad-signal-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    % load CONN data files
    for k=1:length(pathes)
        % load experimental signals
        errorCount = 0;
        subjectNum = 1;
        while errorCount < 2
            sbjfile = sprintf('ROI_Subject%03d_Session001.mat', subjectNum);
            expfile = [base pathes{k} '/conn_project01/data/' sbjfile];

            % check file existance
            if ~exist(expfile, 'file')
                errorCount = errorCount + 1;
                continue;
            end
            % load conn ROI signal file
            load(expfile);
            subjectNum = subjectNum + 1;

            seqLen = size(data{1,START},1);
            si = zeros(ROINUM,seqLen);

            for i=1:ROINUM
                si(i,:) = data{1,START+(i-1)}.';
                roiNames2{1,i} = names{1,START+(i-1)};
            end
            % roi relocation
            si2 = si(ReLOCATE,:);
            signals{end+1} = si2;
            roiNames = roiNames2(1, ReLOCATE);
            
            idx = idx + 1;

            % show signals
            figure;
            plot(si.');
            title([group ':' num2str(ROINUM) ' ROI signals (' num2str(idx) ')']);
        end
    end
    save(outfName, 'signals', 'roiNames');
end

function [weights, meanWeights, stdWeights] = calculateConnectivity(signals, roiNames, group, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    LAG = 3;

    outfName = ['data/ad-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        weights = zeros(ROINUM, ROINUM, length(signals));
        for i=1:length(signals)
            switch(algorithm)
            case 'fc'
                mat = calcFunctionalConnectivity(signals{i});
            case 'gc'
                mat = calcMultivariateGCI(signals{i}, LAG);
            case 'te'
                mat = calcLinueTE(signals{i}, LAG);
            case 'dlcm'
                dlcmName = ['data/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(dlcmName, 'file')
                    load(dlcmName);
                else
                    [si, sig, m, maxsi, minsi] = convert2SigmoidSignal(signals{i});
                    sigLen = size(si,2);
                    inSignal = rand(ROINUM, sigLen);
                    inControl = eye(ROINUM);
                    netDLCM = initDlcmNetwork(si, inSignal, inControl); 
                    % training DLCM network
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
                %            'Plots','training-progress');

                    disp('start training');
                    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
                    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
                    disp(['end training : rsme=' num2str(rsme)]);
                    % calc dlcm-gc
                    mat = calcDlcmGCI(si, inSignal, inControl, netDLCM);
                    
                    save(dlcmName, 'netDLCM', 'si', 'inSignal', 'inControl', 'mat', 'sig', 'm', 'maxsi', 'minsi');
                end
            end
            weights(:,:,i) = mat;
        end
        save(outfName, 'weights', 'roiNames');
    end
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);
    % counting by source region and target region
    meanSource = nanmean(meanWeights,1);
    meanTarget = nanmean(meanWeights,2);
    save(outfName, 'weights', 'meanWeights', 'stdWeights', 'roiNames', 'meanSource', 'meanTarget');

    % show functional conectivity
    figure; 
    switch(algorithm)
    case 'fc'
        clims = [-1,1];
        titleStr = [group ' : Functional Connectivity'];
        sigWeights = meanWeights;
    case 'gc'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : multivariate Granger Causality Index'];
    case 'te'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : Transfer Entropy (LINER)'];
    case 'dlcm'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM Granger Causality Index'];
    end
    imagesc(sigWeights,clims);
    daspect([1 1 1]);
    title(titleStr);
    colorbar;
end

function [normalities, normalitiesP] = calculateNormalityTest(weights, roiNames, group, algorithm)
    % constant value
    ROINUM = size(weights,1);

    outfName = ['data/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-normality.mat'];
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
    % show normality test result
    figure; 
    clims = [0,1];
    imagesc(normalities,clims);
    daspect([1 1 1]);
    title([group '-' algorithm ' normality test result']);
    colorbar;
    % normality test p values
    figure; 
    clims = [0,0.5];
    imagesc(normalitiesP,clims);
    daspect([1 1 1]);
    title([group '-' algorithm ' normality test p values']);
    colorbar;
end


function [utestH, utestP, utestP2] = calculateWilcoxonTest(control, target, roiNames, controlGroup, targetGroup, algorithm)
    % constant value
    ROINUM = size(control,1);

    outfName = ['data/ad-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROINUM) '-utest.mat'];
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
                [p, h] = signrank(x2,y2);
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

    % U test result
    figure; 
    clims = [0,1];
    imagesc(utestH,clims);
    daspect([1 1 1]);
    title([controlGroup '-' targetGroup ' : ' algorithm ' : u test result']);
    colorbar;
    % U test p values
    figure; 
    clims = [0,1];
    imagesc(utestP, clims);
    daspect([1 1 1]);
    title([controlGroup '-' targetGroup ' : ' algorithm ' : u test p values']);
    colorbar;
end
