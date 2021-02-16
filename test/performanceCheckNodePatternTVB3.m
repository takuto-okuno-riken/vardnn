function performanceCheckNodePatternTVB3
    node_nums = [11,22,33,44,55,66];
    num_scan = 55;
    if num_scan == 55  % oh's mouse 98 node. density 0.15. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [1, 1, 1, 1, 1, 1];
    end
    % test sparse and full density
    hz = 64;
    N = 8;

    for i=1:length(node_nums)
        checkingPattern(node_nums(i), num_scan, hz, Gths(i), N, i);
    end
end

function checkingPattern(node_num, num_scan, hz, Gth, N, idx)
    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 12;

    maxLag = 5;
    fname = ['results/rww3/tvb-wongwang3-' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(idx) '-' num2str(hz) 'hz-result.mat'];
    if exist(fname, 'file')
        load(fname);
    else
        % init
        gcAUC = zeros(maxLag,N);
        mvarecAUC = zeros(maxLag,N);
        mpcvarecAUC = zeros(maxLag,N);
        gc2AUC = zeros(maxLag,N);
        mvarec2AUC = zeros(maxLag,N);
        mpcvarec2AUC = zeros(maxLag,N);
        dlAUC = zeros(maxLag,N);
        dlwAUC = zeros(maxLag,N);
        dl2AUC = zeros(maxLag,N);
        dlw2AUC = zeros(maxLag,N);
        dl3AUC = zeros(maxLag,N);
        dlw3AUC = zeros(maxLag,N);
        dl4AUC = zeros(maxLag,N);
        dlw4AUC = zeros(maxLag,N);
        for lags=1:maxLag
            gcROC{lags} = cell(N,2);
            mvarecROC{lags} = cell(N,2);
            mpcvarecROC{lags} = cell(N,2);
            gc2ROC{lags} = cell(N,2);
            mvarec2ROC{lags} = cell(N,2);
            mpcvarec2ROC{lags} = cell(N,2);
            dlROC{lags} = cell(N,2);
            dlwROC{lags} = cell(N,2);
            dl2ROC{lags} = cell(N,2);
            dlw2ROC{lags} = cell(N,2);
            dl3ROC{lags} = cell(N,2);
            dlw3ROC{lags} = cell(N,2);
            dl4ROC{lags} = cell(N,2);
            dlw4ROC{lags} = cell(N,2);
            gcRf{lags} = figure;
            mvarecRf{lags} = figure;
            mpcvarecRf{lags} = figure;
            gc2Rf{lags} = figure;
            mvarec2Rf{lags} = figure;
            mpcvarec2Rf{lags} = figure;
            dlRf{lags} = figure;
            dlwRf{lags} = figure;
            dl2Rf{lags} = figure;
            dlw2Rf{lags} = figure;
            dl3Rf{lags} = figure;
            dlw3Rf{lags} = figure;
            dl4Rf{lags} = figure;
            dlw4Rf{lags} = figure;
        end

        origf = figure;
        origSigf = figure;

        if NumProcessors > 1
            try
                disp('Destroing any existance matlab pool session');
                parpool('close');
            catch
                disp('No matlab pool session found');
            end
            parpool(NumProcessors);
        end

        for k=1:N
            tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(idx) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
            load(tvbFile);
            nodeNum = size(si,1);
            sigLen = size(si,2);

            % show original connection
            figure(origf); plotEC(weights, 'Ground Truth', 1);
            figure(origSigf); plot(t, si);

            % check saved data
            dlcmFile = ['results/rww3/net-patrww3-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(idx) '-' num2str(k) '.mat'];
            netDLCM = cell(maxLag,1);
            netDLCM2 = cell(maxLag,1);
            netDLCM3 = cell(maxLag,1);
            netDLCM4 = cell(maxLag,1);
            if exist(dlcmFile, 'file')
                load(dlcmFile);
            end
            % show DCM signals
            [Y, sig, c, maxsi, minsi] = convert2SigmoidSignal(si);
            [exSignal, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            exControl = eye(nodeNum, nodeNum);
            figure; plot(Y.');
            %figure; plot(exSignal.');

            % training option
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

            % train DLCM with lags
            for lags=1:maxLag
                if isempty(netDLCM{lags})
                    % train DLCM with normal activation function (ReLU)
                    netDLCM{lags} = initDlcmNetwork(Y, exSignal, [], exControl, lags);
                    netDLCM{lags} = trainDlcmNetwork(Y, exSignal, [], exControl, netDLCM{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end
                if isempty(netDLCM2{lags})
                    % train DLCM without activation function (ReLU) (linear case)
                    netDLCM2{lags} = initDlcmNetwork(Y, exSignal, [], exControl, lags, []);
                    netDLCM2{lags} = trainDlcmNetwork(Y, exSignal, [], exControl, netDLCM2{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end
                if isempty(netDLCM3{lags})
                    % train DLCM with normal activation function (ReLU) without exogenous signals
                    netDLCM3{lags} = initDlcmNetwork(Y, [], [], [], lags);
                    netDLCM3{lags} = trainDlcmNetwork(Y, [], [], [], netDLCM3{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end
                if isempty(netDLCM4{lags})
                    % train DLCM without activation function (ReLU) (linear case) without exogenous signals
                    netDLCM4{lags} = initDlcmNetwork(Y, [], [], [], lags, []);
                    netDLCM4{lags} = trainDlcmNetwork(Y, [], [], [], netDLCM4{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end

                % show result of granger causality index (mvGC)
                gcI = calcMultivariateGCI_(si, exSignal, [], exControl, lags);
                figure(gcRf{lags}); hold on; [gcROC{lags}{k,1}, gcROC{lags}{k,2}, gcAUC(lags,k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
                title(['mvGC(' num2str(lags) ')']);

                % show result of granger causality index (mvGC) without exogenous signals
                gcI = calcMultivariateGCI_(si, [], [], [], lags);
                figure(gc2Rf{lags}); hold on; [gc2ROC{lags}{k,1}, gc2ROC{lags}{k,2}, gc2AUC(lags,k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
                title(['mvGC(' num2str(lags) ') (without exogenous)']);

                % extra tests (multivaliate Vector Auto-Regression EC)
                netMVAR = initMvarNetwork(si, exSignal, [], exControl, lags);
                mvarEC = calcMvarEC(netMVAR, [], exControl);
                figure(mvarecRf{lags}); hold on; [mvarecROC{lags}{k,1}, mvarecROC{lags}{k,2}, mvarecAUC(lags,k)] = plotROCcurve(mvarEC, weights, 100, 1, Gth); hold off;
                title(['mVAR(' num2str(lags) ')-EC']);

                % extra tests (multivaliate Vector Auto-Regression EC) without exogenous signals
                netMVAR = initMvarNetwork(si, [], [], [], lags);
                mvarEC = calcMvarEC(netMVAR, [], []);
                figure(mvarec2Rf{lags}); hold on; [mvarec2ROC{lags}{k,1}, mvarec2ROC{lags}{k,2}, mvarec2AUC(lags,k)] = plotROCcurve(mvarEC, weights, 100, 1, Gth); hold off;
                title(['mVAR(' num2str(lags) ')-EC (without exogenous)']);

                % extra tests (multivaliate Principal Component Vector Auto-Regression EC)
                netMVAR = initMpcvarNetwork(si, exSignal, [], exControl, lags);
                mvarEC = calcMpcvarEC(netMVAR, [], exControl);
                figure(mpcvarecRf{lags}); hold on; [mpcvarecROC{lags}{k,1}, mpcvarecROC{lags}{k,2}, mpcvarecAUC(lags,k)] = plotROCcurve(mvarEC, weights, 100, 1, Gth); hold off;
                title(['mPCVAR(' num2str(lags) ')-EC']);

                % extra tests (multivaliate Principal Component Vector Auto-Regression EC) without exogenous signals
                netMVAR = initMpcvarNetwork(si, [], [], [], lags);
                mvarEC = calcMpcvarEC(netMVAR, [], []);
                figure(mpcvarec2Rf{lags}); hold on; [mpcvarec2ROC{lags}{k,1}, mpcvarec2ROC{lags}{k,2}, mpcvarec2AUC(lags,k)] = plotROCcurve(mvarEC, weights, 100, 1, Gth); hold off;
                title(['mPCVAR(' num2str(lags) ')-EC (without exogenous)']);

                % show result of DLCM-GC
                dlGC = calcDlcmGCI(Y, exSignal, [], exControl, netDLCM{lags}, 0);
                figure(dlRf{lags}); hold on; [dlROC{lags}{k,1}, dlROC{lags}{k,2}, dlAUC(lags,k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
                title(['DLCM(' num2str(lags) ')-GC']);

                % show result of DLCM-EC
                dlwEC = calcDlcmEC(netDLCM{lags}, [], exControl, 0);
                figure(dlwRf{lags}); hold on; [dlwROC{lags}{k,1}, dlwROC{lags}{k,2}, dlwAUC(lags,k)] = plotROCcurve(dlwEC, weights, 100, 1, Gth); hold off;
                title(['DLCM(' num2str(lags) ')-EC']);

                % show result of linear DLCM-GC
                dl2GC = calcDlcmGCI(Y, exSignal, [], exControl, netDLCM2{lags}, 0);
                figure(dl2Rf{lags}); hold on; [dl2ROC{lags}{k,1}, dl2ROC{lags}{k,2}, dl2AUC(lags,k)] = plotROCcurve(dl2GC, weights, 100, 1, Gth); hold off;
                title(['linear DLCM(' num2str(lags) ')-GC']);

                % show result of linear DLCM EC 
                dlw2EC = calcDlcmEC(netDLCM2{lags}, [], exControl, 0);
                figure(dlw2Rf{lags}); hold on; [dlw2ROC{lags}{k,1}, dlw2ROC{lags}{k,2}, dlw2AUC(lags,k)] = plotROCcurve(dlw2EC, weights, 100, 1, Gth); hold off;
                title(['linear DLCM(' num2str(lags) ')-EC']);

                % show result of DLCM-GC (without exogenous)
                dl3GC = calcDlcmGCI(Y, [], [], [], netDLCM3{lags}, 0);
                figure(dl3Rf{lags}); hold on; [dl3ROC{lags}{k,1}, dl3ROC{lags}{k,2}, dl3AUC(lags,k)] = plotROCcurve(dl3GC, weights, 100, 1, Gth); hold off;
                title(['DLCM(' num2str(lags) ')-GC (without exogenous)']);

                % show result of DLCM EC (without exogenous)
                dlw3EC = calcDlcmEC(netDLCM3{lags}, [], [], 0);
                figure(dlw3Rf{lags}); hold on; [dlw3ROC{lags}{k,1}, dlw3ROC{lags}{k,2}, dlw3AUC(lags,k)] = plotROCcurve(dlw3EC, weights, 100, 1, Gth); hold off;
                title(['DLCM(' num2str(lags) ')-EC (without exogenous)']);
                
                % show result of linear DLCM-GC (without exogenous)
                dl4GC = calcDlcmGCI(Y, [], [], [], netDLCM4{lags}, 0);
                figure(dl4Rf{lags}); hold on; [dl4ROC{lags}{k,1}, dl4ROC{lags}{k,2}, dl4AUC(lags,k)] = plotROCcurve(dl4GC, weights, 100, 1, Gth); hold off;
                title(['linear DLCM(' num2str(lags) ')-GC (without exogenous)']);

                % show result of linear DLCM EC (without exogenous)
                dlw4EC = calcDlcmEC(netDLCM4{lags}, [], [], 0);
                figure(dlw4Rf{lags}); hold on; [dlw4ROC{lags}{k,1}, dlw4ROC{lags}{k,2}, dlw4AUC(lags,k)] = plotROCcurve(dlw4EC, weights, 100, 1, Gth); hold off;
                title(['linear DLCM(' num2str(lags) ')-EC (without exogenous)']);
            end
        end
        % save result
        save(fname, 'gcAUC',  'gcROC', 'mvarecAUC',  'mvarecROC', 'gc2AUC',  'gc2ROC', 'mvarec2AUC',  'mvarec2ROC', ...
            'mpcvarecAUC', 'mpcvarecROC', 'mpcvarec2AUC', 'mpcvarec2ROC', ...
            'dlAUC', 'dlwAUC', 'dlROC', 'dlwROC', 'dl2AUC', 'dlw2AUC', 'dl2ROC', 'dlw2ROC', ...
            'dl3AUC', 'dlw3AUC', 'dl3ROC', 'dlw3ROC', 'dl4AUC', 'dlw4AUC', 'dl4ROC', 'dlw4ROC');

        % shutdown parallel processing
        if NumProcessors > 1
            delete(gcp('nocreate'))
        end
    end

    % show box plot
    AUCs = nan(N,70);
    r = [1:5];
    AUCs(:,r) = gcAUC.'; r=r+5;
    AUCs(:,r) = gc2AUC.'; r=r+5; % no exogenous
    AUCs(:,r) = mvarecAUC.'; r=r+5;
    AUCs(:,r) = mvarec2AUC.'; r=r+5; % no exogenous
    AUCs(:,r) = mpcvarecAUC.'; r=r+5;
    AUCs(:,r) = mpcvarec2AUC.'; r=r+5; % no exogenous
    AUCs(:,r) = dl2AUC.'; r=r+5;
    AUCs(:,r) = dl4AUC.'; r=r+5;
    AUCs(:,r) = dlw2AUC.'; r=r+5;
    AUCs(:,r) = dlw4AUC.'; r=r+5;
    AUCs(:,r) = dlAUC.'; r=r+5;
    AUCs(:,r) = dl3AUC.'; r=r+5;
    AUCs(:,r) = dlwAUC.'; r=r+5;
    AUCs(:,r) = dlw3AUC.'; r=r+5;
    figure; boxplot(AUCs);
    title(['AUC box plot idx' num2str(idx)]);

    % show average ROC curve of DCM
    figure; 
    hold on;
    for lags=1:maxLag
        plotErrorROCcurve(dlROC{lags}, N, [0.2,0.2,0.2]+(lags*0.1));
        plotErrorROCcurve(dlwROC{lags}, N, [0.2,0.2,0.2]+(lags*0.1));
        plotErrorROCcurve(dl2ROC{lags}, N, [0.2,0.2,0.3]+(lags*0.1));
        plotErrorROCcurve(dlw2ROC{lags}, N, [0.2,0.2,0.3]+(lags*0.1));
        plotAverageROCcurve(dlROC{lags}, N, '--', [0.2,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(dlwROC{lags}, N, '-', [0.2,0.2,0.2]+(lags*0.1),1.0);
        plotAverageROCcurve(dl2ROC{lags}, N, '--', [0.2,0.2,0.3]+(lags*0.1),0.4);
        plotAverageROCcurve(dlw2ROC{lags}, N, '-', [0.2,0.2,0.3]+(lags*0.1),0.4);
    end
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(idx)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
