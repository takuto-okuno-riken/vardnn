% Before using this function, download SPM12 codes from
% https://www.fil.ion.ucl.ac.uk/spm/software/download/
% and add a path "spm12" and sub folders, then remove "spm12/external" folder and sub folders.

% Before using this function, download Dlingam-1.2 codes from
% https://sites.google.com/site/sshimizu06/Dlingamcode
% and add a path "Dlingam-1.2" and sub folders. And also download kernel-ICA 1.2 code from
% https://www.di.ens.fr/~fbach/kernel-ica/index.htm
% and add a path "kernel-ica1_2" and sub folders.

function performanceCheckNodePatternDCM3
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
    disp('network density 0.2'); % self connection
    pP.A = eye(n,n) * 0.2;
    pP.A(5,1) = 0.2 + rand() * 0.3;
    pP.A(8,3) = 0.2 + rand() * 0.3;
    pP.A(6,3) = 0.2 + rand() * 0.3;
    checkingPattern(pP,M,U,N,T,n,TR,options,1);
%%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('network density 0.25');
    pP.A = eye(n,n) * 0.2;
    pP.A(5,1) = 0.3 + rand() * 0.3;
    pP.A(8,3) = 0.3 + rand() * 0.3;
    pP.A(6,3) = 0.3 + rand() * 0.3;
    pP.A(4,6) = 0.3 + rand() * 0.3;
    pP.A(7,6) = 0.3 + rand() * 0.3;
    pP.A(4,8) = 0.3 + rand() * 0.3;
    checkingPattern(pP,M,U,N,T,n,TR,options,2);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('network density 0.3');
    pP.A = eye(n,n) * 0.2;
    pP.A = addPattern6(pP.A,0.3,0.2);
    checkingPattern(pP,M,U,N,T,n,TR,options,6);
%%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('network density 0.41');
    pP.A = eye(n,n) * 0.15;
    pP.A = addPattern6(pP.A,0.3,0.2);
    pP.A = addPattern7(pP.A);
    checkingPattern(pP,M,U,N,T,n,TR,options,7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('network density 0.5');
    pP.A = eye(n,n) * 0.1;
    pP.A = addPattern6(pP.A,0.2,0.2);
    pP.A = addPattern7(pP.A);
    pP.A = addPattern8(pP.A);
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
function A = addPattern7(A)
    A(1,4) = 0.1 + rand() * 0.2;
    A(3,4) = 0.1 + rand() * 0.2;
    A(5,4) = 0.1 + rand() * 0.2;
    A(6,4) = 0.1 + rand() * 0.2;
    A(3,5) = 0.1 + rand() * 0.2;
    A(6,5) = 0.1 + rand() * 0.2;
end
function A = addPattern8(A)
    A(2,4) = 0.1 + rand() * 0.2;
    A(2,5) = 0.1 + rand() * 0.2;
    A(2,6) = 0.1 + rand() * 0.2;
    A(2,7) = 0.1 + rand() * 0.2;
    A(2,8) = 0.1 + rand() * 0.2;
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
function [FC, dlGC, gcI] = checkingPattern(pP,M,U,N,T,n,TR,options,idx)
    % show original connection
    figure; plotDcmEC(pP.A);

    fcAUC = zeros(1,N);
    pcAUC = zeros(1,N);
    pcpcAUC = zeros(1,N);
    lsopcAUC = zeros(1,N);
    plspcAUC = zeros(1,N);
    wcsAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    dlmAUC = zeros(1,N);
    dlgAUC = zeros(1,N);
    pcdlAUC = zeros(1,N);
    pcdlwAUC = zeros(1,N);
    svlpcAUC = zeros(1,N);
    svgpcAUC = zeros(1,N);
    svrpcAUC = zeros(1,N);
    gppcAUC = zeros(1,N);
    trpcAUC = zeros(1,N);
    rfpcAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    pcpcROC = cell(N,2);
    lsopcROC = cell(N,2);
    plspcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);
    dlmROC = cell(N,2);
    dlgROC = cell(N,2);
    pcdlROC = cell(N,2);
    pcdlwROC = cell(N,2);
    svlpcROC = cell(N,2);
    svgpcROC = cell(N,2);
    svrpcROC = cell(N,2);
    gppcROC = cell(N,2);
    trpcROC = cell(N,2);
    rfpcROC = cell(N,2);
    fcRf = figure;
    pcRf = figure;
    pcpcRf = figure;
    lsopcRf = figure;
    plspcRf = figure;
    wcsRf = figure;
    gcRf = figure;
    pgcRf = figure;
    dlRf = figure;
    dlwRf = figure;
    dlmRf = figure;
    dlgRf = figure;
    pcdlRf = figure;
    pcdlwRf = figure;
    svlpcRf = figure;
    svgpcRf = figure;
    svrpcRf = figure;
    gppcRf = figure;
    trpcRf = figure;
    rfpcRf = figure;

    fname = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    load(fname);

    % calc input signal and node BOLD signals
    for k=1:N
        netFile = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        netDLCM = [];
        if exist(netFile, 'file')
            load(netFile);
        else
            % generate signal by DCM
            U.u = spm_rand_mar(T+50,n,1/2)/8;       % endogenous fluctuations
            y   = spm_int_J(pP,M,U);                % integrate with observer
            y2  = y(51:end,:);
            u2  = U.u(51:end,:);
            si = y2.';
            data = si; % for RNN-GC
            save(netFile, 'netDLCM', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data');
        end

        % show result of FC
        fg = figure; FC = plotFunctionalConnectivity(y2.'); close(fg);
        figure(fcRf); hold on; [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = plotROCcurve(FC, pP.A); hold off;
        title('FC');
        % show result of PC
        fg = figure; PC = plotPartialCorrelation(y2.'); close(fg);
        figure(pcRf); hold on; [pcROC{k,1}, pcROC{k,2}, pcAUC(k)] = plotROCcurve(PC, pP.A); hold off;
        title('PC');
        % PC and PLSPC diff check
        if isempty(plspcROC{k,1})
            PC2 = calcPLSPartialCorrelation(y2.'); % calc PLS PC
            Z = PC - PC2; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-PLSPC=' num2str(pcdiff)]);
    %        figure; clims = [-1 1]; imagesc(Z,clims); title('PC - PLSPC');
            figure(plspcRf); hold on; [plspcROC{k,1}, plspcROC{k,2}, plspcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('PLS-PC');
        end
        % show result of PCA-PC
        if isempty(pcpcROC{k,1})
            PC2 = calcPcPartialCorrelation(y2.'); % calc PCA PC
            Z = PC - PC2; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-PCAPC=' num2str(pcdiff)]);
    %        figure; clims = [-1 1]; imagesc(Z,clims); title('PC - PCAPC');
            figure(pcpcRf); hold on; [pcpcROC{k,1}, pcpcROC{k,2}, pcpcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('PCA-PC');
        end
        % show result of Lasso-PC
        if isempty(lsopcROC{k,1})
            [lambda, alpha, errMat] = estimateLassoParamsForPC(y2.', [], [], [], 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
            PC2 = calcLassoPartialCorrelation(y2.', [], [], [], lambda, alpha); % calc Lasso PC
            Z = PC - PC2; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-LassoPC=' num2str(pcdiff)]);
    %        figure; clims = [-1 1]; imagesc(Z,clims); title('PC - LassoPC');
            figure(lsopcRf); hold on; [lsopcROC{k,1}, lsopcROC{k,2}, lsopcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('Lasso-PC');
        end
        % show result of SVR-PC
        if isempty(svlpcROC{k,1})
            PC2 = calcSvPartialCorrelation(y2.',[], [], [], 'linear');
            figure(svlpcRf); hold on; [svlpcROC{k,1}, svlpcROC{k,2}, svlpcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('SVl-PC');
        end
        if isempty(svgpcROC{k,1})
            PC2 = calcSvPartialCorrelation(y2.',[], [], [], 'gaussian');
            figure(svgpcRf); hold on; [svgpcROC{k,1}, svgpcROC{k,2}, svgpcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('SVg-PC');
        end
        if isempty(svrpcROC{k,1})
            PC2 = calcSvPartialCorrelation(y2.',[], [], [], 'rbf');
            figure(svrpcRf); hold on; [svrpcROC{k,1}, svrpcROC{k,2}, svrpcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('SVr-PC');
        end
        % show result of GP-PC
        if isempty(gppcROC{k,1})
            PC2 = calcGpPartialCorrelation(y2.');
            figure(gppcRf); hold on; [gppcROC{k,1}, gppcROC{k,2}, gppcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('GP-PC');
        end
        % show result of Tree-PC
        if isempty(trpcROC{k,1})
            PC2 = calcTreePartialCorrelation(y2.');
            figure(trpcRf); hold on; [trpcROC{k,1}, trpcROC{k,2}, trpcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('Tree-PC');
        end
        % show result of RF-PC
        if isempty(rfpcROC{k,1})
            PC2 = calcRfPartialCorrelation(y2.');
            figure(rfpcRf); hold on; [rfpcROC{k,1}, rfpcROC{k,2}, rfpcAUC(k)] = plotROCcurve(PC2, pP.A); hold off;
            title('RF-PC');
        end

        % show result of WCS
        if isempty(wcsROC{k,1})
            fg = figure; WCS = plotWaveletCoherence(y2.'); close(fg);
            figure(wcsRf); hold on; [wcsROC{k,1}, wcsROC{k,2}, wcsAUC(k)] = plotROCcurve(WCS, pP.A); hold off;
            title('WCS');
        end
        % show result of granger causality index (mvGC)
        if isempty(gcROC{k,1})
            fg = figure; gcI = plotMultivariateGCI_(y2.',[],[],[],3,0); close(fg);
            figure(gcRf); hold on; [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = plotROCcurve(gcI, pP.A); hold off;
            title('mvGC');
        end
        % mvGC and mPLSVARGC diff check
        %{
        netMPLSVAR = initMplsvarNetwork(y2.', [], [], [], 3);
        mplsvarGC = calcMplsvarGCI(y2.', [], [], [], netMPLSVAR);
        Z = gcI - mplsvarGC; gcdiff=nanmean(abs(Z),'all'); disp(['mae of mvGC-mplsvarGC=' num2str(gcdiff) ' / ' num2str(max(max(gcI)))]);
        figure; clims = [-1 1]; imagesc(Z,clims); title('mvGC - mplsvarGC');
        %}
        % show result of granger causality index (pwGC)
        if isempty(pgcROC{k,1})
            fg = figure; gcI = plotPairwiseGCI(y2.',[],[],[],3,0); close(fg);
            figure(pgcRf); hold on; [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = plotROCcurve(gcI, pP.A); hold off;
            title('pwGC');
        end
        % show result of DirectLiNGAM
        if isempty(dlgROC{k,1})
            fg = figure; Aest = plotDirectLiNGAM(y2.'); close(fg);
            figure(dlgRf); hold on; [dlgROC{k,1}, dlgROC{k,2}, dlgAUC(k)] = plotROCcurve(Aest, pP.A); hold off;
            title('dLiNGAM');
        end
        % show DCM signals
        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(y2.', 0);
        [exSignal, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(u2.', 0);
        % si = si - 0.5; (bad VARDNN-DI)
        % exSignal = exSignal - 0.5; (bad VARDNN-DI)
        % si = y2.'; % test raw data (bad VARDNN-DI)
        % exSignal = u2.'; % test raw data (bad VARDNN-DI)
        % si = y2.' - min(y2,[],'all'); % test raw data (nice VARDNN-DI)
        % exSignal = u2.' - min(u2,[],'all');; % test raw data (nice VARDNN-DI)
        exControl = eye(n,n);
        figure; plot(si.');
        %figure; plot(exSignal.');

        % train VARDNN
        nodeNum = size(si,1);
        sigLen = size(si,2);
        if isempty(netDLCM)
            % layer parameters
            netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
            % training VARDNN network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 3);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.1, ...
                'Verbose',false);
        %            'Plots','training-progress');

            disp('start training');
            netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
            % recoverty training
            %[netDLCM, time] = recoveryTrainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
            save(netFile, 'netDLCM', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
        end

        % show signals after training
        figure; [S, t,mae,maeerr] = plotPredictSignals(si,exSignal,[],exControl,netDLCM);
        disp(['t=' num2str(t) ', mae=' num2str(mae)]);

        % show result of VARDNN-GC
        if isempty(dlROC{k,1})
            fg = figure; dlGC = plotMvarDnnGCI(si, exSignal, [], exControl, netDLCM, 0); close(fg);
            figure(dlRf); hold on; [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = plotROCcurve(dlGC, pP.A); hold off;
            title('VARDNN-GC');
        end
        % show result of VARDNN weight causality index (VARDNN-WCI) as VARDNN-DI
        if isempty(dlwROC{k,1})
            fg = figure; dlw = plotMvarDnnDI(netDLCM, [], exControl, 0); close(fg);
            figure(dlwRf); hold on; [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = plotROCcurve(dlw, pP.A); hold off;
            title('VARDNN-DI');
        end
        % show result of VARDNN mean impact value (VARDNN-MIV)
        if isempty(dlmROC{k,1})
            fg = figure; [~,dlm] = calcMvarDnnMIV(si, exSignal, [], exControl, netDLCM, 0); close(fg);
            figure(dlmRf); hold on; [dlmROC{k,1}, dlmROC{k,2}, dlmAUC(k)] = plotROCcurve(dlm, pP.A); hold off;
            title('VARDNN-MIV');
        end
        
        % train PC-VARDNN
        pcvarFile = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) 'pcvar.mat'];
        netMPC = [];
        if exist(pcvarFile, 'file')
            load(pcvarFile);
        else
            % layer parameters
            netMPC = initMpcvarDnnNetwork(si, exSignal, [], exControl);
            % training PC-VARDNN network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 3);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.1, ...
                'Verbose',false);

            disp('start training');
            netMPC = trainMpcvarDnnNetwork(si, exSignal, [], exControl, netMPC, options);
            save(pcvarFile, 'netMPC');
        end

        % show result of PC-VARDNN-GC
        if isempty(pcdlROC{k,1})
            fg = figure; dlGC = plotMpcvarDnnGCI(si, exSignal, [], exControl, netMPC, 0); close(fg);
            figure(pcdlRf); hold on; [pcdlROC{k,1}, pcdlROC{k,2}, pcdlAUC(k)] = plotROCcurve(dlGC, pP.A); hold off;
            title('PC-VARDNN-GC');
        end
        % show result of PC-VARDNN weight causality index (PC-VARDNN-WCI) as PC-VARDNN-DI
        if isempty(pcdlwROC{k,1})
            fg = figure; dlw = plotMpcvarDnnDI(netMPC, [], exControl, 0); close(fg);
            figure(pcdlwRf); hold on; [pcdlwROC{k,1}, pcdlwROC{k,2}, pcdlwAUC(k)] = plotROCcurve(dlw, pP.A); hold off;
            title('PC-VARDNN-DI');
        end
    end
    
    % save result
    if exist('nvdiAUC','var')
        save(fname, 'fcAUC','pcAUC','pcpcAUC','plspcAUC','lsopcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','dlmAUC','dlgAUC','pcdlAUC','pcdlwAUC','dcmAUC', ...
            'rnnAUC','linueAUC','nnnueAUC','pcsAUC','cpcAUC','fgesAUC','fcaAUC','tsfcAUC','tsfcaAUC', ...
            'mvarecAUC','pvarecAUC','mpcvarecAUC','mpcvargcAUC','ppcvarecAUC','ppcvargcAUC',...
            'mplsecAUC','mplsgcAUC','pplsecAUC','pplsgcAUC',...
            'plsoecAUC','mlsoecAUC','plsogcAUC','mlsogcAUC','pcgcAUC','dls1AUC','dls3AUC','msvmdiAUC','msvmgcAUC','mgpdiAUC','mgpgcAUC','mgpediAUC',...
            'nvdiAUC','nvmiAUC','trdiAUC','trmiAUC','rfdiAUC','rfmiAUC', ...
            'svlpcAUC','svgpcAUC','svrpcAUC','gppcAUC','trpcAUC','rfpcAUC', ...
            'fcROC','pcROC','pcpcROC','plspcROC','lsopcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlmROC','dlgROC','pcdlROC','pcdlwROC','dcmROC', ...
            'rnnROC','linueROC','nnnueROC','pcsROC','cpcROC','fgesROC','fcaROC','tsfcROC','tsfcaROC', ...
            'mvarecROC','pvarecROC','mpcvarecROC','mpcvargcROC','ppcvarecROC','ppcvargcROC', ...
            'mplsecROC','mplsgcROC','pplsecROC','pplsgcROC', ...
            'plsoecROC','mlsoecROC','plsogcROC','mlsogcROC','pcgcROC','dls1ROC','dls3ROC','msvmdiROC','msvmgcROC','mgpdiROC','mgpgcROC','mgpediROC', ...
            'nvdiROC','nvmiROC','trdiROC','trmiROC','rfdiROC','rfmiROC', ...
            'svlpcROC','svgpcROC','svrpcROC','gppcROC','trpcROC','rfpcROC');
    else
        save(fname, 'fcAUC', 'pcAUC', 'pcpcAUC', 'plspcAUC', 'lsopcAUC', 'wcsAUC', 'gcAUC', 'pgcAUC', 'dlAUC', 'dlwAUC', 'dlmAUC', 'dlgAUC', 'pcdlAUC', 'pcdlwAUC', ...
            'fcROC','pcROC','pcpcROC','plspcROC','lsopcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlmROC','dlgROC','pcdlROC','pcdlwROC');
    end
end

