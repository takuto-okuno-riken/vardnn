% this function is only for ADNI2 alzheimer analysis

function [weights, meanWeights, stdWeights] = calculateConnectivity(signals, roiNames, group, algorithm, rawFlag)
    if nargin < 5
        rawFlag = 0;
    end
    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 14;

    % constant value
    ROINUM = size(signals{1},1);
    LAG = 3;

    outfName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        if NumProcessors > 1
            try
                disp('Destroing any existance matlab pool session');
                parpool('close');
            catch
                disp('No matlab pool session found');
            end
            parpool(NumProcessors);
        end

        weights = zeros(ROINUM, ROINUM, length(signals));
        parfor i=1:length(signals)    % for parallel processing
%        for i=1:length(signals)
            switch(algorithm)
            case 'fc'
                mat = calcFunctionalConnectivity(signals{i});
            case 'pc'
                mat = calcPartialCorrelation(signals{i});
            case 'wcs'
                fName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(fName, 'file')
                    f = load(fName);
                    mat = f.mat;
                else
                    mat = calcWaveletCoherence(signals{i});
                    parsavemat(fName, mat);
                end
            case 'gc'
                mat = calcMultivariateGCI(signals{i}, LAG);
            case 'pgc'
                mat = calcPairwiseGCI(signals{i}, LAG);
            case 'te'
                mat = calcLinueTE(signals{i}, LAG);
            case 'pcs'
                csvFile = ['results/tetrad/pcs-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'cpc'
                csvFile = ['results/tetrad/cpc-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'fges'
                csvFile = ['results/tetrad/fges-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'dlg'
                fName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(fName, 'file')
                    f = load(fName);
                    mat = f.mat;
                else
                    mat = calcDirectLiNGAM(signals{i});
                    parsavemat(fName, mat);
                end
            case 'dlcm'
                dlcmName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(dlcmName, 'file')
                    f = load(dlcmName);
                    mat = f.mat;
                else
                    if rawFlag
                        si = signals{i};
                        sig=0; c=0; maxsi=0; minsi=0;
                    else
                        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
                    end
                    % si = signals{i} - nanmin(signals{i}, [], 'all'); % simple linear transform
                    sigLen = size(si,2);
                    inSignal = rand(ROINUM, sigLen);
                    inControl = eye(ROINUM);
                    netDLCM = initDlcmNetwork(si, inSignal, [], inControl); 
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
                    netDLCM = trainDlcmNetwork(si, inSignal, [], inControl, netDLCM, options);
                    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
                    disp(['end training : rsme=' num2str(rsme)]);
                    % calc dlcm-gc
                    mat = calcDlcmGCI(si, inSignal, [], inControl, netDLCM);
                    
                    parsavedlsm(dlcmName, netDLCM, si, inSignal, inControl, mat, sig, c, maxsi, minsi);
                end
            case 'dlcmrc' % should be called after dlcm
                outName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(outName, 'file')
                    f = load(outName);
                    mat = f.mat;
                else
                    dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                    f = load(dlcmName);
                    if isfield(f,'c'), c=f.c; else c=f.m; end
                    % recoverty training
                    options = trainingOptions('adam', ...
                        'ExecutionEnvironment','cpu', ...
                        'MaxEpochs', 1000, ...
                        'MiniBatchSize',ceil(size(f.si,2) / 3), ...
                        'Shuffle','every-epoch', ...
                        'GradientThreshold',5,...
                        'L2Regularization',0.05, ...
                        'Verbose',false);
                    [netDLCM, time] = recoveryTrainDlcmNetwork(f.si, f.inSignal, [], f.inControl, f.netDLCM, options);
                    mat = calcDlcmGCI(f.si, f.inSignal, [], f.inControl, netDLCM);
                    parsavedlsm(outName, netDLCM, f.si, f.inSignal, f.inControl, mat, f.sig, c, f.maxsi, f.minsi);
                end
            case {'dlw','dlwrc'} % should be called after dlcm
                if strcmp(algorithm, 'dlw')
                    dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                else
                    dlcmName = ['results/ad-dlcmrc-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                end
                f = load(dlcmName);
                mat = calcDlcmEC(f.netDLCM, [], f.inControl);
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
    case 'pc'
        clims = [-1,1];
        titleStr = [group ' : Partial Correlation'];
        sigWeights = meanWeights;
    case 'wcs'
        clims = [-1,1];
        titleStr = [group ' : Wavelet Coherence Cross Spectrum'];
        sigWeights = meanWeights;
    case 'gc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : multivariate Granger Causality Index'];
    case 'pgc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pairwise Granger Causality Index'];
    case 'te'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : Transfer Entropy (LINER)'];
    case 'pcs'
        clims = [-1,1];
        titleStr = [group ' : PC-stable-max'];
        sigWeights = meanWeights;
    case 'cpc'
        clims = [-1,1];
        titleStr = [group ' : Conservative PC'];
        sigWeights = meanWeights;
    case 'fges'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : FGES'];
    case 'dlg'
        clims = [-1,1];
        titleStr = [group ' : Direct LiNGAM'];
        sigWeights = meanWeights;
    case {'dlcm','dlcmrc'}
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM Granger Causality Index'];
    case {'dlw','dlwrc'}
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM Weight Causality Index'];
    end
    imagesc(sigWeights,clims);
    daspect([1 1 1]);
    title(titleStr);
    colorbar;

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
end

function parsavemat(fName, mat)
    save(fName, 'mat');
end

function parsavedlsm(dlcmName, netDLCM, si, inSignal, inControl, mat, sig, c, maxsi, minsi)
    save(dlcmName, 'netDLCM', 'si', 'inSignal', 'inControl', 'mat', 'sig', 'c', 'maxsi', 'minsi');
end
