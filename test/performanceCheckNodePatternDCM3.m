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
    wcsAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    dlgAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);
    dlgROC = cell(N,2);
    fcRf = figure;
    pcRf = figure;
    wcsRf = figure;
    gcRf = figure;
    pgcRf = figure;
    dlRf = figure;
    dlwRf = figure;
    dlgRf = figure;

    % calc input signal and node BOLD signals
    for k=1:N
        dlcmFile = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        netDLCM = [];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            % generate signal by DCM
            U.u = spm_rand_mar(T+50,n,1/2)/8;       % endogenous fluctuations
            y   = spm_int_J(pP,M,U);                % integrate with observer
            y2  = y(51:end,:);
            u2  = U.u(51:end,:);
            si = y2.';
            data = si; % for RNN-GC
            save(dlcmFile, 'netDLCM', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data');
        end

        % show result of FC
        fg = figure; FC = plotFunctionalConnectivity(y2.'); close(fg);
        figure(fcRf); hold on; [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = plotROCcurve(FC, pP.A); hold off;
        % show result of PC
        fg = figure; PC = plotPartialCorrelation(y2.'); close(fg);
        figure(pcRf); hold on; [pcROC{k,1}, pcROC{k,2}, pcAUC(k)] = plotROCcurve(PC, pP.A); hold off;
        % show result of WCS
        fg = figure; WCS = plotWaveletCoherence(y2.'); close(fg);
        figure(wcsRf); hold on; [wcsROC{k,1}, wcsROC{k,2}, wcsAUC(k)] = plotROCcurve(WCS, pP.A); hold off;
        % show result of granger causality index (mvGC)
        fg = figure; gcI = plotMultivariateGCI_(y2.',[],[],[],3,0); close(fg);
        figure(gcRf); hold on; [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = plotROCcurve(gcI, pP.A); hold off;
        title('mvGC');
        % show result of granger causality index (pwGC)
        fg = figure; gcI = plotPairwiseGCI(y2.',[],[],[],3,0); close(fg);
        figure(pgcRf); hold on; [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = plotROCcurve(gcI, pP.A); hold off;
        title('pwGC');
        % show result of DirectLiNGAM
        fg = figure; Aest = plotDirectLiNGAM(y2.'); close(fg);
        figure(dlgRf); hold on; [dlgROC{k,1}, dlgROC{k,2}, dlgAUC(k)] = plotROCcurve(Aest, pP.A); hold off;
        title('dLiNGAM');

        % show DCM signals
        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(y2.', 0);
        [exSignal, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(u2.', 0);
        % si = si - 0.5; (bad DLCM-EC)
        % exSignal = exSignal - 0.5; (bad DLCM-EC)
        % si = y2.'; % test raw data (bad DLCM-EC)
        % exSignal = u2.'; % test raw data (bad DLCM-EC)
        % si = y2.' - min(y2,[],'all'); % test raw data (nice DLCM-EC)
        % exSignal = u2.' - min(u2,[],'all');; % test raw data (nice DLCM-EC)
        exControl = eye(n,n);
        figure; plot(si.');
        %figure; plot(exSignal.');

        % train DLCM
        nodeNum = size(si,1);
        sigLen = size(si,2);
        if isempty(netDLCM)
            % layer parameters
            netDLCM = initDlcmNetwork(si, exSignal, [], exControl);
            % training DLCM network
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
            netDLCM = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
            % recoverty training
            %[netDLCM, time] = recoveryTrainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
            save(dlcmFile, 'netDLCM', 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'si', 'data', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
        end

        % show signals after training
        figure; [S, t,mae,maeerr] = plotPredictSignals(si,exSignal,[],exControl,netDLCM);
        disp(['t=' num2str(t) ', mae=' num2str(mae)]);

        % show result of DLCM-GC
        fg = figure; dlGC = plotDlcmGCI(si, exSignal, [], exControl, netDLCM, 0); close(fg);
        figure(dlRf); hold on; [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = plotROCcurve(dlGC, pP.A); hold off;
        title('DLCM-GC');

        % show result of DLCM weight causality index (DLCM-wci) as DLCM-EC
        fg = figure; dlwGC = plotDlcmEC(netDLCM, [], exControl, 0); close(fg);
        figure(dlwRf); hold on; [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = plotROCcurve(dlwGC, pP.A); hold off;
        title('DLCM-EC');
    end
    fname = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    save(fname, 'fcAUC', 'pcAUC', 'wcsAUC', 'gcAUC', 'pgcAUC', 'dlAUC', 'dlwAUC', 'dlgAUC', 'fcROC','pcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlgROC');
end

