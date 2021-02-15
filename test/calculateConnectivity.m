% this function is only for ADNI2 alzheimer analysis

function [weights, meanWeights, stdWeights, subweights] = calculateConnectivity(signals, roiNames, group, algorithm, isRaw, lags, isAutoExo, activateFunc)
    if nargin < 8, activateFunc = @reluLayer; end
    if nargin < 7, isAutoExo = 0; end
    if nargin < 6, lags = 1; end
    if nargin < 5, isRaw = 0; end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 11;

    % constant value
    ROINUM = size(signals{1},1);
    sigLen = size(signals{1},2);

    weights = zeros(ROINUM, ROINUM, length(signals));
    subweights = zeros(ROINUM, ROINUM+1, length(signals));

%    lagpat = ["gc","pgc","te","tsfc","tsfca","mvarec","dlcm","dlw"];
%    if lags>1 && contains(algorithm,lagpat), lagStr=num2str(lags); else lagStr=''; end
    if lags>1, lagStr=num2str(lags); else lagStr=''; end
    if isAutoExo>0, exoStr='_ex'; else exoStr=''; end
    if isempty(activateFunc), linStr='_lin'; else linStr=''; end
    outfName = ['results/ad-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '.mat'];
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

%        parfor i=1:length(signals)    % for parallel processing
        for i=1:length(signals)
            if isAutoExo > 0
                exSignal = rand(ROINUM, sigLen);
                exControl = eye(ROINUM);
            else
                exSignal = [];
                exControl = [];
            end
            switch(algorithm)
            case 'fc'
                mat = calcFunctionalConnectivity(signals{i});
            case 'fca'
                mat = calcFunctionalConnectivityAbs(signals{i});
            case 'tsfc'
                mat = calcTimeShiftedCorrelation(signals{i}, exSignal, [], exControl, lags);
            case 'tsfca'
                mat = calcTimeShiftedCorrelationAbs(signals{i}, exSignal, [], exControl, lags);
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
                mat = calcMultivariateGCI_(signals{i}, exSignal, [], exControl, lags);
                %mat = calcMultivariateGCI(signals{i}, exSignal, [], exControl, lags);
            case 'pgc'
                mat = calcPairwiseGCI(signals{i}, exSignal, [], exControl, lags);
            case 'te'
                mat = calcLinueTE(signals{i}, exSignal, [], exControl, lags);
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
            case 'mvarec'
                netMVAR = initMvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcMvarEC(netMVAR, [], exControl); % |Zi-Zi\j| version
            case 'mvar'
                netMVAR = initMvarNetwork(signals{i}, exSignal, [], exControl, lags);
                [mat, sub] = calcMvarEC(netMVAR, [], exControl); % |Zi-Zi\j| version
                mat = repmat(sub(:,1), [1 size(mat,2)]) - sub(:,2:end); % subtract Zi-Zi\j version
            case 'mpcvarec'
                netMPCVAR = initMpcvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcMpcvarEC(netMPCVAR, [], exControl);
            case 'dlcm'
                dlcmName = ['results/ad-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(dlcmName, 'file')
                    f = load(dlcmName);
                    mat = f.mat;
                else
                    if isRaw
                        si = signals{i};
                        sig=0; c=0; maxsi=0; minsi=0;
                    else
                        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
                    end
                    % si = signals{i} - nanmin(signals{i}, [], 'all'); % simple linear transform
                    netDLCM = initDlcmNetwork(si, exSignal, [], exControl, lags, activateFunc); 
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
                    netDLCM = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
                    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
                    disp(['end training : rsme=' num2str(rsme)]);
                    % calc dlcm-gc
                    mat = calcDlcmGCI(si, exSignal, [], exControl, netDLCM);
                    
                    parsavedlsm(dlcmName, netDLCM, si, exSignal, exControl, mat, sig, c, maxsi, minsi);
                end
            case 'dlcmrc' % should be called after dlcm
                outName = ['results/ad-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(outName, 'file')
                    f = load(outName);
                    mat = f.mat;
                else
                    dlcmName = ['results/ad-dlcm' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                    f = load(dlcmName);
                    if isfield(f,'c'), c=f.c; else c=f.m; end % for compatibility
                    if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
                    if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
                    % recoverty training
                    options = trainingOptions('adam', ...
                        'ExecutionEnvironment','cpu', ...
                        'MaxEpochs', 1000, ...
                        'MiniBatchSize',ceil(size(f.si,2) / 3), ...
                        'Shuffle','every-epoch', ...
                        'GradientThreshold',5,...
                        'L2Regularization',0.05, ...
                        'Verbose',false);
                    [netDLCM, time] = recoveryTrainDlcmNetwork(f.si, f.exSignal, [], f.exControl, f.netDLCM, options);
                    mat = calcDlcmGCI(f.si, f.exSignal, [], f.exControl, netDLCM);
                    parsavedlsm(outName, netDLCM, f.si, f.exSignal, f.exControl, mat, f.sig, c, f.maxsi, f.minsi);
                end
            case {'dlw','dlwrc'} % should be called after dlcm
                if strcmp(algorithm, 'dlw')
                    dlcmName = ['results/ad-dlcm' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                else
                    dlcmName = ['results/ad-dlcmrc' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                end
                f = load(dlcmName);
                if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
                [mat, subweights(:,:,i)] = calcDlcmEC(f.netDLCM, [], f.exControl);
            end
            weights(:,:,i) = mat;
        end
        save(outfName, 'weights', 'roiNames', 'subweights');
    end
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);
    % counting by source region and target region
    meanSource = nanmean(meanWeights,1);
    meanTarget = nanmean(meanWeights,2);
    save(outfName, 'weights', 'meanWeights', 'stdWeights', 'roiNames', 'meanSource', 'meanTarget', 'subweights');

    % show functional conectivity
    figure; 
    switch(algorithm)
    case 'fc'
        clims = [-1,1];
        titleStr = [group ' : Functional Connectivity'];
        sigWeights = meanWeights;
    case 'fca'
        clims = [0,1];
        titleStr = [group ' : Functional Connectivity (Abs)'];
        sigWeights = meanWeights;
    case 'tsfc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : Time shifted Functional Connectivity (' num2str(lags) ')'];
    case 'tsfca'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : Time shifted Functional Connectivity (Abs) (' num2str(lags) ')'];
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
        titleStr = [group ' : multivariate Granger Causality(' num2str(lags) ') Index'];
    case 'pgc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pairwise Granger Causality(' num2str(lags) ') Index'];
    case 'te'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : Transfer Entropy (LINER) (' num2str(lags) ')'];
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
        titleStr = [group ' : DLCM(' num2str(lags) ') Granger Causality Index'];
    case {'dlw','dlwrc'}
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM(' num2str(lags) ') Weight Causality Index'];
    case 'mvarec'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : mVAR(' num2str(lags) ')-EC Index'];
    case 'mvar'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : mVAR(' num2str(lags) ') Index'];
    case 'mpcvarec'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : mPCVAR(' num2str(lags) ')-EC Index'];
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

function parsavedlsm(dlcmName, netDLCM, si, exSignal, exControl, mat, sig, c, maxsi, minsi)
    save(dlcmName, 'netDLCM', 'si', 'exSignal', 'exControl', 'mat', 'sig', 'c', 'maxsi', 'minsi');
end
