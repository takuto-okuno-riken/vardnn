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
    n  = 16;                               % number of regions or nodes
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
    LAG = 3;   % lag time for GC

    % initialize parameters to save file or load previous one
    FCtime = zeros(maxK,N);
    GCtime = zeros(maxK,N);
    DLtime = zeros(maxK,N);

    % load spectrum DCM performance check result
    fname = ['results/DCM_demo-rand' num2str(n) '-' num2str(N) 'x' num2str(maxK) '.mat'];
    load(fname);

    for k = 1:maxK
        % get BOLD signal of DCM simulation result
        pP  = pPs{k};     % A,B,C,D params of minimum error row
        U.u = Uus{k,N};  % endogenous fluctuations
        y   = spm_int_J(pP,M,U);            % integrate with observer
        si  = y(ySt:end,:).';
        inSignal = U.u(ySt:end,:).';
        inControl = eye(n,n);

        data = si;
        dlcmFile = ['results/net-timeD-' num2str(n) '-' num2str(N) 'x' num2str(k) '.mat'];
        save(dlcmFile, 'data','inSignal');

        % calculate FC and get computational time
        ticHdl = tic;
        mat = calcFunctionalConnectivity(si);
        time = toc(ticHdl);
        FCtime(k,1) = time;
        disp(['finish FC! t = ' num2str(time) 's']);

        % calculate GC and get computational time
        ticHdl = tic;
        mat = calcMultivariateGCI(si, LAG);
        time = toc(ticHdl);
        GCtime(k,1) = time;
        disp(['finish GC! t = ' num2str(time) 's']);

        % normalize signal to [0, 1] 
        si = convert2SigmoidSignal(si,0);
        inSignal = convert2SigmoidSignal(inSignal,0);

        % do training or load DLCM network
        dlcmFile = ['results/net-time-' num2str(n) '-' num2str(N) 'x' num2str(k) '.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            ticHdl = tic;
            % init DLCM network
            netDLCM = initDlcmNetwork(si, inSignal, [], inControl);

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

            % training DLCM network
            netDLCM = trainDlcmNetwork(si, inSignal, [], inControl, netDLCM, options);
            % calc dlcm-gc
            mat = calcDlcmGCI(si, inSignal, [], inControl, netDLCM);
            time = toc(ticHdl);
            disp(['finish calculating DLCM-GC! t = ' num2str(time) 's']);

            DLtime(k,1) = time;

            %plotDlcmWeight(netDLCM);
            save(dlcmFile, 'netDLCM','mat','DLtime');
        end
    end
    % save result
%    dlcmFile = ['results/net-time-' num2str(n) '-' num2str(N) 'x' num2str(k) '-result.mat'];
%    save(dlcmFile, 'DLtime', 'FCtime', 'GCtime');
end


