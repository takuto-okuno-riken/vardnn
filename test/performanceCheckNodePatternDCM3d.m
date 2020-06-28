
% this script should run after performanceCheckNodePatternDCM3
function performanceCheckNodePatternDCM3d
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
    disp('network density 0.191');
    checkingPattern(pP,M,U,N,T,n,TR,options,1);
%%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('network density 0.25');
    checkingPattern(pP,M,U,N,T,n,TR,options,2);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('network density 0.304');
    checkingPattern(pP,M,U,N,T,n,TR,options,6);
%%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('network density 0.411');
    checkingPattern(pP,M,U,N,T,n,TR,options,7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('network density 0.5');
    checkingPattern(pP,M,U,N,T,n,TR,options,8);
%%}
end

%% 
function checkingPattern(pP,M,U,N,T,n,TR,options,idx)
    fname = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    load(fname);

    % init
    dcmAUC = zeros(1,N);
    dcmROC = cell(N,2);
    dcmRf = figure;
    origf = figure;

    % calc input signal and node BOLD signals
    for k=1:N
        dlcmFile = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(dlcmFile);

        % show original connection
        figure(origf); plotDcmEC(pP.A);

        dcmFile = ['results/net-pat3-'  num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) 'd.mat'];
        if exist(dcmFile, 'file')
            load(dcmFile);
        else
            CSD = {};
            RMS = [];
        end
        orgUu = U.u;

        % initialize DCM stcuct
        DCM = struct();
        DCM.options = options;

        DCM.a    = ones(n,n);
        DCM.b    = zeros(n,n,0);
        DCM.c    = eye(n,n);
        DCM.d    = zeros(n,n,0);

        DCM.Y.dt = TR;
        DCM.U.dt = TR;

        % performance check of DCM inversion
        for j=length(CSD)+1:6
            % generate signal by DCM
            % add white noise to get different BOLD signal. Only DCM does get several pattern 
            % of input signal for estimating one connectivity to get nicer ROC curve :(
            U.u = orgUu + spm_rand_mar(T+50,n,1/2)/64; % endogenous fluctuations
            y   = spm_int_J(pP,M,U);                % integrate with observer
            y2  = y(51:end,:);
            u2  = U.u(51:end,:);

            % response
            % -----------------------------------------------------------------
            DCM.Y.y  = y2;
            DCM.U.u  = u2;

            % nonlinear system identification (Variational Laplace)
            % =================================================================
            CSD{end + 1} = spm_dcm_fmri_csd(DCM);
            BPA          = spm_dcm_average(CSD,'simulation',1);

            dp   = BPA.Ep.A - pP.A;
            dp   = dp - diag(diag(dp));
            RMS(end + 1) = sqrt(mean(dp(~~dp).^2))

            save(dcmFile, 'pP', 'M', 'U','n','TR', 'y2', 'u2', 'RMS', 'CSD');
        end

        % show estimated A by DCM
        BPA = spm_dcm_average(CSD,'simulation',1);
        figure; plotDcmEC(BPA.Ep.A,0);
        figure(dcmRf); hold on; [dcmROC{k,1}, dcmROC{k,2}, dcmAUC(k)] = plotROCcurve(BPA.Ep.A, pP.A); hold off;
    end
    save(fname, 'fcAUC','gcAUC','dlAUC','dcmAUC', 'fcROC','gcROC','dlROC','dcmROC');
end

