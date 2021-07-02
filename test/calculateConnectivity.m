% this function is only for ADNI2 alzheimer analysis

function [weights, meanWeights, stdWeights, subweights] = calculateConnectivity(signals, roiNames, group, algorithm, isRaw, lags, isAutoExo, activateFunc, sigmoidAlpha)
    if nargin < 9, sigmoidAlpha = 1; end
    if nargin < 8, activateFunc = @reluLayer; end
    if nargin < 7, isAutoExo = 0; end
    if nargin < 6, lags = 1; end
    if nargin < 5, isRaw = 0; end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 11;

    % constant value
    ROINUM = size(signals{1},1);

    weights = zeros(ROINUM, ROINUM, length(signals));
    subweights = zeros(ROINUM, ROINUM+1, length(signals));

    global resultsPath;
    global resultsPrefix;
%    lagpat = ["gc","pgc","te","tsfc","tsfca","mvarec","dlcm","dlw"];
%    if lags>1 && contains(algorithm,lagpat), lagStr=num2str(lags); else lagStr=''; end
    if lags>1, lagStr=num2str(lags); else lagStr=''; end
    if isAutoExo>0, exoStr='_ex'; else exoStr=''; end
    if isempty(activateFunc), linStr='_lin'; else linStr=''; end
    outfName = [resultsPath '/' resultsPrefix '-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '.mat'];
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

        parfor i=1:length(signals)    % for parallel processing
%        for i=1:length(signals)
            sigLen = size(signals{i},2);
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
%                mat = calcPartialCorrelation(signals{i});
                mat = calcPartialCorrelation__(signals{i});
            case 'plspc'
                mat = calcPLSPartialCorrelation(signals{i});
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
                csvFile = [resultsPath '/tetrad/pcs-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'cpc'
                csvFile = [resultsPath '/tetrad/cpc-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'fges'
                csvFile = [resultsPath '/tetrad/fges-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'dlg'
                fName = [resultsPath '/' resultsPrefix '-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
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
                [~, ~, mat] = calcMvarEC(netMVAR, [], exControl); % |Zi-Zi\j| version
            case 'pvarec'
                mat = calcPvarEC(signals{i}, exSignal, [], exControl, lags); % |Zi-Zi\j| version
            case 'mpcvarec'
                netMPCVAR = initMpcvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcMpcvarEC(netMPCVAR, [], exControl);
            case 'mpcvargc'
                netMPCVAR = initMpcvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcMpcvarGCI(signals{i}, exSignal, [], exControl, netMPCVAR);
            case 'ppcvarec'
                netPPCVAR = initPpcvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcPpcvarEC(netPPCVAR, [], exControl);
            case 'ppcvargc'
                netPPCVAR = initPpcvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcPpcvarGCI(signals{i}, exSignal, [], exControl, netPPCVAR);
            case 'mplsvarec'
                netMPLSVAR = initMplsvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcMplsvarEC(netMPLSVAR, [], exControl);
            case 'mplsvargc'
                netMPLSVAR = initMplsvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcMplsvarGCI(signals{i}, exSignal, [], exControl, netMPLSVAR);
            case 'pplsvarec'
                netPPLSVAR = initPplsvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcPplsvarEC(netPPLSVAR, [], exControl);
            case 'pplsvargc'
                netPPLSVAR = initPplsvarNetwork(signals{i}, exSignal, [], exControl, lags);
                mat = calcPplsvarGCI(signals{i}, exSignal, [], exControl, netPPLSVAR);
            case {'dlcm', 'dlcmB'}
                dlcmName = [resultsPath '/' resultsPrefix '-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(dlcmName, 'file')
                    f = load(dlcmName);
                    mat = f.mat;
                else
                    if isRaw
                        si = signals{i};
                        sig=0; c=0; maxsi=0; minsi=0;
                    else
                        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i},NaN,sigmoidAlpha);
                    end
                    % si = signals{i} - nanmin(signals{i}, [], 'all'); % simple linear transform
                    netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl, lags, activateFunc); 
                    % training VARDNN network
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
                    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
                    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
                    disp(['end training : rsme=' num2str(rsme)]);
                    % calc dlcm-gc
                    mat = calcMvarDnnGCI(si, exSignal, [], exControl, netDLCM);
                    
                    parsavedlsm(dlcmName, netDLCM, si, exSignal, exControl, mat, sig, c, maxsi, minsi);
                end
            case 'dlcmrc' % should be called after dlcm
                outName = [resultsPath '/' resultsPrefix '-' algorithm lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(outName, 'file')
                    f = load(outName);
                    mat = f.mat;
                else
                    dlcmName = [resultsPath '/' resultsPrefix '-dlcm' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
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
                    mat = calcMvarDnnGCI(f.si, f.exSignal, [], f.exControl, netDLCM);
                    parsavedlsm(outName, netDLCM, f.si, f.exSignal, f.exControl, mat, f.sig, c, f.maxsi, f.minsi);
                end
            case {'dlw','dlwB','dlwrc'} % should be called after dlcm
                if strcmp(algorithm, 'dlw')
                    dlcmName = [resultsPath '/' resultsPrefix '-dlcm' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                elseif strcmp(algorithm, 'dlwB')
                    dlcmName = [resultsPath '/' resultsPrefix '-dlcmB' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                else
                    dlcmName = [resultsPath '/' resultsPrefix '-dlcmrc' lagStr exoStr linStr '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                end
                f = load(dlcmName);
                if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
                [mat, subweights(:,:,i)] = calcMvarDnnEC(f.netDLCM, [], f.exControl);
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
    case 'plspc'
        clims = [-1,1];
        titleStr = [group ' : PLS Partial Correlation'];
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
    case {'dlcm','dlcmB','dlcmrc'}
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM(' num2str(lags) ') Granger Causality Index'];
    case {'dlw','dlwB','dlwrc'}
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
    case 'pvarec'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pVAR(' num2str(lags) ')-EC Index'];
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
    case 'mpcvargc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : mPCVAR(' num2str(lags) ')-GC Index'];
    case 'ppcvarec'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pPCVAR(' num2str(lags) ')-EC Index'];
    case 'ppcvargc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pPCVAR(' num2str(lags) ')-GC Index'];
    case 'mplsvarec'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : mPLSVAR(' num2str(lags) ')-EC Index'];
    case 'mplsvargc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : mPLSVAR(' num2str(lags) ')-GC Index'];
    case 'pplsvarec'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pPLSVAR(' num2str(lags) ')-EC Index'];
    case 'pplsvargc'
        sigma = std(meanWeights(:),1,'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pPLSVAR(' num2str(lags) ')-GC Index'];
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
