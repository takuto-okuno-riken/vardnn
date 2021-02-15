% Before using this function, download SPM12 codes from
% https://www.fil.ion.ucl.ac.uk/spm/software/download/
% and add a path "spm12" and sub folders, then remove "spm12/external" folder and sub folders.

function performanceCheckNodePatternDCM6
    % set global random stream and shuffle it
    myStream=RandStream('mt19937ar');
    RandStream.setGlobalStream(myStream);
    rng('shuffle');

    % DEM Structure: create random inputs
    % -------------------------------------------------------------------------
    N  = 8;
    T  = 300;                             % number of observations (scans)
    TR = 2;                               % repetition time or timing
    n  = 8;                               % number of regions or nodes
    t  = (1:T)*TR;                        % observation times

    % priors
    % -------------------------------------------------------------------------
    options.maxnodes   = 4;  % effective number of nodes, 4 is better than n

    options.nonlinear  = 0;
    options.two_state  = 0;
    options.stochastic = 0;
    options.centre     = 1;
    options.induced    = 1;

    A   = eye(n,n);
    B   = zeros(n,n,0);
    C   = zeros(n,n);
    D   = zeros(n,n,0);
    pP  = spm_dcm_fmri_priors(A,B,C,D,options);

    pP.C = eye(n,n);
    pP.transit = randn(n,1)/16;

    % integrate states
    % -------------------------------------------------------------------------
    U.dt = TR;
    M.f  = 'spm_fx_fmri';
    M.x  = sparse(n,5);
    M.g   = 'spm_gx_fmri';

    %% pattern 1 -------------------------------------------------
%%{
    disp('network density 0.05');
    pP.A = rand(n,n)/5 - 0.1;
    pP.A(5,1) = 0.5 + rand() * 0.2;
    pP.A(8,3) = 0.5 + rand() * 0.2;
    pP.A(6,3) = 0.5 + rand() * 0.2;
    checkingPattern(pP,M,U,N,T,n,TR,options,1);
%%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('network density 0.11');
    pP.A = rand(n,n)/5 - 0.1;
    pP.A(5,1) = 0.5 + rand() * 0.2;
    pP.A(8,3) = 0.5 + rand() * 0.2;
    pP.A(6,3) = 0.5 + rand() * 0.2;
    pP.A(4,6) = 0.5 + rand() * 0.2;
    pP.A(7,6) = 0.5 + rand() * 0.2;
    pP.A(4,8) = 0.5 + rand() * 0.2;
    checkingPattern(pP,M,U,N,T,n,TR,options,2);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('network density 0.16');
    pP.A = rand(n,n)/5 - 0.1;
    pP.A = addPattern6(pP.A,0.3,0.2);
    checkingPattern(pP,M,U,N,T,n,TR,options,6);
%%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('network density 0.27');
    pP.A = rand(n,n)/5 - 0.1;
    pP.A = addPattern6(pP.A,0.25,0.1);
    pP.A = addPattern7(pP.A,0.25,0.2);
    checkingPattern(pP,M,U,N,T,n,TR,options,7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('network density 0.36');
    pP.A = rand(n,n)/5 - 0.1;
    pP.A = addPattern6(pP.A,0.25,0.1);
    pP.A = addPattern7(pP.A,0.25,0.2);
    pP.A = addPattern8(pP.A,0.25,0.2);
    checkingPattern(pP,M,U,N,T,n,TR,options,8);
%%}
end

function A = addPattern6(A,base,range)
    A(3,1) = base + rand() * range;
    A(5,1) = base + rand() * range;
    A(8,3) = base + rand() * range;
    A(8,5) = base + rand() * range;
    A(5,8) = base + rand() * range;
    A(7,4) = base + rand() * range;
    A(6,7) = base + rand() * range;
    A(7,6) = base + rand() * range;
    A(6,8) = base + rand() * range;
end
function A = addPattern7(A,base,range)
    A(1,4) = base + rand() * range;
    A(3,4) = base + rand() * range;
    A(5,4) = base + rand() * range;
    A(6,4) = base + rand() * range;
    A(3,5) = base + rand() * range;
    A(6,5) = base + rand() * range;
end
function A = addPattern8(A,base,range)
    A(2,4) = base + rand() * range;
    A(2,5) = base + rand() * range;
    A(2,6) = base + rand() * range;
    A(2,7) = base + rand() * range;
    A(2,8) = base + rand() * range;
end
function A = addPatternD(A,n)
    A(2,:) = NaN;
    A(:,2) = NaN;
    didx = find(eye(n,n)>0);
    A(didx) = NaN;
    idx = find(A>0);
    len = length(idx);
    idx = idx(randperm(len));
    A(idx(1:len-20)) = 0;
    idx = find(isnan(A));
    A(idx) = 0;
    A = A + eye(n,n) * (0.2);
end

%% 
function [dlGC] = checkingPattern(pP,M,U,N,T,n,TR,options,idx)
    % show original connection
    figure; plotDcmEC(pP.A);
    maxLag = 5;

    fname = ['results/net-pat6-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    if exist(fname, 'file')
        load(fname);
    else
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

        % calc input signal and node BOLD signals
        for k=1:N
            % read same y2, u2 si signals
            pat4File = ['results/net-pat4-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
            load(pat4File);

            % check saved data
            dlcmFile = ['results/net-pat6-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
            netDLCM = cell(maxLag,1);
            netDLCM2 = cell(maxLag,1);
            netDLCM3 = cell(maxLag,1);
            netDLCM4 = cell(maxLag,1);
            if exist(dlcmFile, 'file')
                load(dlcmFile);
            end
            % show DCM signals
            [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(y2.', 0);
            [exSignal, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(u2.', 0);
            exControl = eye(n,n);
            figure; plot(si.');
            %figure; plot(exSignal.');

            % training option
            sigLen = size(si,2);
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
                    netDLCM{lags} = initDlcmNetwork(si, exSignal, [], exControl, lags);
                    netDLCM{lags} = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end
                if isempty(netDLCM2{lags})
                    % train DLCM without activation function (ReLU) (linear case)
                    netDLCM2{lags} = initDlcmNetwork(si, exSignal, [], exControl, lags, []);
                    netDLCM2{lags} = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM2{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end
                if isempty(netDLCM3{lags})
                    % train DLCM with normal activation function (ReLU) without exogenous signals
                    netDLCM3{lags} = initDlcmNetwork(si, [], [], [], lags);
                    netDLCM3{lags} = trainDlcmNetwork(si, [], [], [], netDLCM3{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end
                if isempty(netDLCM4{lags})
                    % train DLCM without activation function (ReLU) (linear case) without exogenous signals
                    netDLCM4{lags} = initDlcmNetwork(si, [], [], [], lags, []);
                    netDLCM4{lags} = trainDlcmNetwork(si, [], [], [], netDLCM4{lags}, options);
                    save(dlcmFile, 'netDLCM', 'netDLCM2', 'netDLCM3', 'netDLCM4', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
                end

                % show result of granger causality index (mvGC)
                gcI = calcMultivariateGCI_(y2.', exSignal, [], exControl, lags);
                figure(gcRf{lags}); hold on; [gcROC{lags}{k,1}, gcROC{lags}{k,2}, gcAUC(lags,k)] = plotROCcurve(gcI, pP.A, 100, 1, 0.2); hold off;
                title(['mvGC(' num2str(lags) ')']);

                % show result of granger causality index (mvGC) without exogenous signals
                gcI = calcMultivariateGCI_(y2.', [], [], [], lags);
                figure(gc2Rf{lags}); hold on; [gc2ROC{lags}{k,1}, gc2ROC{lags}{k,2}, gc2AUC(lags,k)] = plotROCcurve(gcI, pP.A, 100, 1, 0.2); hold off;
                title(['mvGC(' num2str(lags) ') (without exogenous)']);

                % extra tests (multivaliate Vector Auto-Regression EC)
                netMVAR = initMvarNetwork(y2.', exSignal, [], exControl, lags);
                mvarEC = calcMvarEC(netMVAR, [], exControl);
                figure(mvarecRf{lags}); hold on; [mvarecROC{lags}{k,1}, mvarecROC{lags}{k,2}, mvarecAUC(lags,k)] = plotROCcurve(mvarEC, pP.A, 100, 1, 0.2); hold off;
                title(['mVAR(' num2str(lags) ')-EC']);

                % extra tests (multivaliate Vector Auto-Regression EC) without exogenous signals
                netMVAR = initMvarNetwork(y2.', [], [], [], lags);
                mvarEC = calcMvarEC(netMVAR, [], []);
                figure(mvarec2Rf{lags}); hold on; [mvarec2ROC{lags}{k,1}, mvarec2ROC{lags}{k,2}, mvarec2AUC(lags,k)] = plotROCcurve(mvarEC, pP.A, 100, 1, 0.2); hold off;
                title(['mVAR(' num2str(lags) ')-EC (without exogenous)']);

                % extra tests (multivaliate Principal Component Vector Auto-Regression EC)
                netMVAR = initMpcvarNetwork(y2.', exSignal, [], exControl, lags);
                mvarEC = calcMpcvarEC(netMVAR, [], exControl);
                figure(mpcvarecRf{lags}); hold on; [mpcvarecROC{lags}{k,1}, mpcvarecROC{lags}{k,2}, mpcvarecAUC(lags,k)] = plotROCcurve(mvarEC, pP.A, 100, 1, 0.2); hold off;
                title(['mPCVAR(' num2str(lags) ')-EC']);

                % extra tests (multivaliate Principal Component Vector Auto-Regression EC) without exogenous signals
                netMVAR = initMpcvarNetwork(y2.', [], [], [], lags);
                mvarEC = calcMpcvarEC(netMVAR, [], []);
                figure(mpcvarec2Rf{lags}); hold on; [mpcvarec2ROC{lags}{k,1}, mpcvarec2ROC{lags}{k,2}, mpcvarec2AUC(lags,k)] = plotROCcurve(mvarEC, pP.A, 100, 1, 0.2); hold off;
                title(['mPCVAR(' num2str(lags) ')-EC (without exogenous)']);

                % show result of DLCM-GC
                dlGC = calcDlcmGCI(si, exSignal, [], exControl, netDLCM{lags}, 0);
                figure(dlRf{lags}); hold on; [dlROC{lags}{k,1}, dlROC{lags}{k,2}, dlAUC(lags,k)] = plotROCcurve(dlGC, pP.A, 100, 1, 0.2); hold off;
                title(['DLCM(' num2str(lags) ')-GC']);

                % show result of DLCM-EC
                dlwEC = calcDlcmEC(netDLCM{lags}, [], exControl, 0);
                figure(dlwRf{lags}); hold on; [dlwROC{lags}{k,1}, dlwROC{lags}{k,2}, dlwAUC(lags,k)] = plotROCcurve(dlwEC, pP.A, 100, 1, 0.2); hold off;
                title(['DLCM(' num2str(lags) ')-EC']);

                % show result of linear DLCM-GC
                dl2GC = calcDlcmGCI(si, exSignal, [], exControl, netDLCM2{lags}, 0);
                figure(dl2Rf{lags}); hold on; [dl2ROC{lags}{k,1}, dl2ROC{lags}{k,2}, dl2AUC(lags,k)] = plotROCcurve(dl2GC, pP.A, 100, 1, 0.2); hold off;
                title(['linear DLCM(' num2str(lags) ')-GC']);

                % show result of linear DLCM EC 
                dlw2EC = calcDlcmEC(netDLCM2{lags}, [], exControl, 0);
                figure(dlw2Rf{lags}); hold on; [dlw2ROC{lags}{k,1}, dlw2ROC{lags}{k,2}, dlw2AUC(lags,k)] = plotROCcurve(dlw2EC, pP.A, 100, 1, 0.2); hold off;
                title(['linear DLCM(' num2str(lags) ')-EC']);

                % show result of DLCM-GC (without exogenous)
                dl3GC = calcDlcmGCI(si, [], [], [], netDLCM3{lags}, 0);
                figure(dl3Rf{lags}); hold on; [dl3ROC{lags}{k,1}, dl3ROC{lags}{k,2}, dl3AUC(lags,k)] = plotROCcurve(dl3GC, pP.A, 100, 1, 0.2); hold off;
                title(['DLCM(' num2str(lags) ')-GC (without exogenous)']);

                % show result of DLCM EC (without exogenous)
                dlw3EC = calcDlcmEC(netDLCM3{lags}, [], [], 0);
                figure(dlw3Rf{lags}); hold on; [dlw3ROC{lags}{k,1}, dlw3ROC{lags}{k,2}, dlw3AUC(lags,k)] = plotROCcurve(dlw3EC, pP.A, 100, 1, 0.2); hold off;
                title(['DLCM(' num2str(lags) ')-EC (without exogenous)']);
                
                % show result of linear DLCM-GC (without exogenous)
                dl4GC = calcDlcmGCI(si, [], [], [], netDLCM4{lags}, 0);
                figure(dl4Rf{lags}); hold on; [dl4ROC{lags}{k,1}, dl4ROC{lags}{k,2}, dl4AUC(lags,k)] = plotROCcurve(dl4GC, pP.A, 100, 1, 0.2); hold off;
                title(['linear DLCM(' num2str(lags) ')-GC (without exogenous)']);

                % show result of linear DLCM EC (without exogenous)
                dlw4EC = calcDlcmEC(netDLCM4{lags}, [], [], 0);
                figure(dlw4Rf{lags}); hold on; [dlw4ROC{lags}{k,1}, dlw4ROC{lags}{k,2}, dlw4AUC(lags,k)] = plotROCcurve(dlw4EC, pP.A, 100, 1, 0.2); hold off;
                title(['linear DLCM(' num2str(lags) ')-EC (without exogenous)']);
            end
        end
        save(fname, 'gcAUC',  'gcROC', 'mvarecAUC',  'mvarecROC', 'gc2AUC',  'gc2ROC', 'mvarec2AUC',  'mvarec2ROC', ...
            'mpcvarecAUC', 'mpcvarecROC', 'mpcvarec2AUC', 'mpcvarec2ROC', ...
            'dlAUC', 'dlwAUC', 'dlROC', 'dlwROC', 'dl2AUC', 'dlw2AUC', 'dl2ROC', 'dlw2ROC', ...
            'dl3AUC', 'dlw3AUC', 'dl3ROC', 'dlw3ROC', 'dl4AUC', 'dlw4AUC', 'dl4ROC', 'dlw4ROC');
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

