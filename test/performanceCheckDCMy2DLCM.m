%==========================================================================
% train and simulate DLCM by simulated BOLD signal of DCM
% Before running this script, please add path of spm12 directory
%
function performanceCheckDCMy2DLCM

    % DEM Structure: create random inputs
    % -------------------------------------------------------------------------
    N  = 12;                              % number of runs
    T  = 200;                             % number of observations (scans)
    TR = 2;                               % repetition time or timing
    n  = 12;                               % number of regions or nodes
    t  = (1:T)*TR;                        % observation times
    
    % integrate states
    % -------------------------------------------------------------------------
    U.dt = TR;
    M.f  = 'spm_fx_fmri';
    M.x  = sparse(n,5);
    M.g   = 'spm_gx_fmri';

    ySt = 30; % BOLD signal starting point (for DCM inversion)
    T = 199 + ySt;

    maxK = 10; % number of trial

    % initialize parameters to save file or load previous one
    DLMAEs = zeros(maxK,N);
    DLinvTms = zeros(maxK,N);
    DLCorr = zeros(maxK,N);

    % load spectrum DCM performance check result
    fname = ['results/DCM_demo-rand' num2str(n) '-' num2str(N) 'x' num2str(maxK) '.mat'];
    load(fname);

    for k = 1:maxK
        % get BOLD signal of DCM simulation result
        pP  = pPs{k};     % A,B,C,D params of minimum error row
        U.u = Uus{k,N};  % endogenous fluctuations
        y   = spm_int_J(pP,M,U);            % integrate with observer
        si  = y(ySt:end,:).';
        exSignal = U.u(ySt:end,:).';
        exControl = eye(n,n);

        % normalize signal to [0, 1] 
        si = bold2dnnSignal(si);
        exSignal = bold2dnnSignal(exSignal);

        % do training or load DLCM network
        dlcmFile = ['results/net-vsDCM-' num2str(n) '-' num2str(N) 'x' num2str(k) '.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            % init VARDNN network
            netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);

            % set training options
            maxEpochs = 1000;
            sigLen = size(si,2);
            miniBatchSize = ceil(sigLen / 3);

            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'L2Regularization',0.1, ...
                'GradientThreshold',5,...
                'Verbose',false);
        %            'Plots','training-progress');
%                'GradientThresholdMethod', 'global-l2norm' , ...

            % training VARDNN network
            netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
            % recover training 
            [netDLCM, time, mae] = recoveryTrainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
            disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

            DLMAEs(k,1) = mae;
            DLinvTms(k,1) = netDLCM.trainTime + netDLCM.recoverTrainTime;

            %plotMvarDnnWeight(netDLCM);
            save(dlcmFile, 'netDLCM', 'DLMAEs', 'DLinvTms', 'DLCorr');
        end
        
        % simulate DLCM network with 1st frame & exogenous input signal
        [S, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, netDLCM);

        % show original & simulated signal correlation        
        figure; DLCorr(k,1) = plotTwoSignalsCorrelation(S, si);
        save(dlcmFile, 'netDLCM', 'DLMAEs', 'DLinvTms', 'DLCorr');
    end
end


