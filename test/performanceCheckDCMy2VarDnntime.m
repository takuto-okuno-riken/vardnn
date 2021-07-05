%==========================================================================
% train and simulate VARDNN by simulated BOLD signal of DCM
% Before running this script, please add path of spm12 directory
%
function performanceCheckDCMy2VarDnntime

    % DEM Structure: create random inputs
    % -------------------------------------------------------------------------
    N  = 12;                              % number of runs
    T  = 200;                             % number of observations (scans)
    TR = 2;                               % repetition time or timing
    n  = 8;                               % number of regions or nodes
    t  = (1:T)*TR;                        % observation times
    
    % integrate states
    % -------------------------------------------------------------------------
    U.dt = TR;
    M.f  = 'spm_fx_fmri';
    M.x  = sparse(n,5);
    M.g   = 'spm_gx_fmri';

    ySt = 30; % BOLD signal starting point (for DCM inversion)
    T = 199 + ySt;

    maxK = 20; % number of trial
    LAG = 3;   % lag time for GC

    % initialize parameters to save file or load previous one
    FCtime = zeros(maxK,N);
    PCtime = zeros(maxK,N);
    WCtime = zeros(maxK,N);
    mGCtime = zeros(maxK,N);
    pGCtime = zeros(maxK,N);
    DLtime = zeros(maxK,N);
    DLEtime = zeros(maxK,N);
    dLiNtime = zeros(maxK,N);

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

        data = si;
        netFile = ['results/net-timeD-' num2str(n) '-' num2str(N) 'x' num2str(k) '.mat'];
        save(netFile, 'data','exSignal');

        % calculate FC and get computational time
        ticHdl = tic;
        mat = calcFunctionalConnectivity(si);
        time = toc(ticHdl);
        FCtime(k,1) = time;
        disp(['finish FC! t = ' num2str(time) 's']);

        % calculate PC and get computational time
        ticHdl = tic;
        mat = calcPartialCorrelation(si);
        time = toc(ticHdl);
        PCtime(k,1) = time;
        disp(['finish PC! t = ' num2str(time) 's']);

        % calculate WC and get computational time
        ticHdl = tic;
        mat = calcWaveletCoherence(si);
        time = toc(ticHdl);
        WCtime(k,1) = time;
        disp(['finish WC! t = ' num2str(time) 's']);

        % calculate mvGC and get computational time
        ticHdl = tic;
        mat = calcMultivariateGCI_(si, exSignal, [], exControl, LAG);
        time = toc(ticHdl);
        mGCtime(k,1) = time;
        disp(['finish mvGC! t = ' num2str(time) 's']);

        % calculate pwGC and get computational time
        ticHdl = tic;
        mat = calcPairwiseGCI(si, exSignal, [], exControl, LAG);
        time = toc(ticHdl);
        pGCtime(k,1) = time;
        disp(['finish pwGC! t = ' num2str(time) 's']);

        % calculate dLiN and get computational time
        ticHdl = tic;
        mat = calcDirectLiNGAM(si);
        time = toc(ticHdl);
        dLiNtime(k,1) = time;
        disp(['finish dLiNGAM! t = ' num2str(time) 's']);
        
        % normalize signal to [0, 1] 
        si = convert2SigmoidSignal(si,0);
        exSignal = convert2SigmoidSignal(exSignal,0);

        % do training or load VARDNN network
        netFile = ['results/net-time-' num2str(n) '-' num2str(N) 'x' num2str(k) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            ticHdl = tic;
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
            % calc VARDNN-GC
            mat = calcMvarDnnGCI(si, exSignal, [], exControl, netDLCM);
            time = toc(ticHdl);
            disp(['finish calculating VARDNN-GC! t = ' num2str(time) 's']);

            DLtime(k,1) = time;

            %plotMvarDnnWeight(netDLCM);
            save(netFile, 'netDLCM','mat','DLtime');
        end
        
        % calculate VARDNN-DI and get computational time
        ticHdl = tic;
        mat = calcMvarDnnEC(netDLCM, [], exControl);
        time = toc(ticHdl) + netDLCM.trainTime;
        DLEtime(k,1) = time;
        disp(['finish VARDNN-DI! t = ' num2str(time) 's']);        
    end

    % save result
    timeFile = ['results/net-time-' num2str(n) '-' num2str(N) 'x' num2str(k) '-result.mat'];
    save(timeFile, 'FCtime', 'PCtime', 'WCtime', 'mGCtime', 'pGCtime', 'DLtime', 'DLEtime', 'dLiNtime');
end


