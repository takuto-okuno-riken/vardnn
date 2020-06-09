%==========================================================================
% performance check spectrum DCM
% Before running this script, please add path of spm12 directory
% This script is modified by 'DEM_demo_large_fMRI' function that is
% Demonstration of DCM for CSD (fMRI) with simulated responses in spm12
%
function performanceCheckSpm12DCMfMRI
    rng('default')

    % DEM Structure: create random inputs
    % -------------------------------------------------------------------------
    N  = 12;                               % number of runs
    T  = 200;                             % number of observations (scans)
    TR = 2;                               % repetition time or timing
    n  = 8;                               % number of regions or nodes
    t  = (1:T)*TR;                        % observation times

    % priors
    % -------------------------------------------------------------------------
    options.maxnodes   = 4;               % effective number of nodes

    options.nonlinear  = 0;
    options.two_state  = 0;
    options.stochastic = 0;
    options.centre     = 1;
    options.induced    = 1;

    A   = ones(n,n);
    B   = zeros(n,n,0);
    C   = zeros(n,n);
    D   = zeros(n,n,0);
    pP  = spm_dcm_fmri_priors(A,B,C,D,options);

    pP.C = eye(n,n);
    pP.transit = randn(n,1)/16;

    % simulate response to endogenous fluctuations
    %==========================================================================

    % integrate states
    % -------------------------------------------------------------------------
    U.dt = TR;
    M.f  = 'spm_fx_fmri';
    M.x  = sparse(n,5);
    M.g   = 'spm_gx_fmri';

    maxK = 20; % number of trial

    % initialize parameters to save file or load previous one
    fname = ['performance_check/DCM_demo-rand' num2str(n) '-' num2str(N) 'x' num2str(maxK) '.mat'];
    if exist(fname, 'file')
        load(fname);
        Corr = zeros(maxK,N);
    else
        pPs = cell(maxK,1);
        Uus = cell(maxK,N);
        BPAs = cell(maxK,N);
        RMSs = zeros(maxK,N);
        MAEs = zeros(maxK,N);
        invTms = zeros(maxK,N);
        Corr = zeros(maxK,N);
    end
    % find restart point
    idx = find(invTms.'==0);
    if isempty(idx)
        startK = maxK + 1;
    else
        startK = floor((idx(1)-1)/N + 1);
    end
startK = 1;
    % performance check of DCM inversion
    for k=startK:maxK

        ySt = 30; % BOLD signal starting point (for DCM inversion)
        T = 199 + ySt;

        % initialize DCM stcuct
        DCM = struct();
        DCM.options = options;

        DCM.a    = ones(n,n);
        DCM.b    = zeros(n,n,0);
        DCM.c    = zeros(n,1);
        DCM.d    = zeros(n,n,0);

        DCM.Y.dt = TR;
        DCM.U.u  = zeros(T,1);
        DCM.U.dt = TR;

        % find starting point
        rowTms = invTms(k,:);
        startI = find(rowTms==0);
startI(1) = 1;
        % generate original connectivty with posterior expectations and generating
        % new data for Bayesian parameter averaging
        %==========================================================================
        pP.A = rand(n,n)/4 - 0.125;
        orgpP = pP;
        pPs{k} = pP;

        if startI(1)==1, CSD   = {}; end
        RMS   = [];
        Qp    = [];
        Pp    = pP.A - diag(diag(pP.A));

        for i = startI(1):N
            try
                % integrate states
                % -----------------------------------------------------------------
                if isempty(Uus{k,i})
                    U.u = spm_rand_mar(T,n,1/2)/8;      % endogenous fluctuations
                    Uus{k,i} = U.u; % keep input & params to save
                else
                    U.u = Uus{k,i};
                end
                y    = spm_int_J(pP,M,U);            % integrate with observer
                orgY = y;

                % response
                % -----------------------------------------------------------------
                DCM.Y.y  = y(ySt:end,:);
%                dcmTicH = tic;
                % nonlinear system identification (Variational Laplace)
                % =================================================================
                CSD{end + 1} = spm_dcm_fmri_csd(DCM);
                BPA          = spm_dcm_average(CSD,'simulation',1);
%                invTms(k,i) = toc(dcmTicH);

                BPAs{k,i} = BPA; % keep input & params to save
%{
                % MAP estimates
                % -----------------------------------------------------------------
                qp   = CSD{end}.Ep.A;
                Qp(:,end + 1) = spm_vec(qp - diag(diag(qp)));

                % root mean square error
                % -----------------------------------------------------------------
                dp   = BPA.Ep.A - pP.A;
                dp   = dp - diag(diag(dp));
                RMS(end + 1) = sqrt(mean(dp(~~dp).^2));

                RMSs(k,i) = RMS(end);

                % summary
                % -----------------------------------------------------------------
                spm_figure('Getwin','Figure 2'); clf

                subplot(2,1,1); hold off
                spm_plot_ci(BPA.Ep.A(:),BPA.Cp(1:n*n,1:n*n)), hold on
                bar(pP.A(:),1/4), hold off
                title('True and MAP connections','FontSize',16)
                axis square

                subplot(2,2,3); cla
                plot(Pp(:),Qp,'cd','MarkerSize',8),hold on
                plot(Pp,BPA.Ep.A - diag(diag(BPA.Ep.A)),'b.','MarkerSize',16), hold off
                title('True and MAP connections (Extrinsic)','FontSize',16)
                xlabel('True')
                ylabel('Estimated')
                axis square

                subplot(2,2,4);
                plot(RMS)
                title('root mean square error','FontSize',16)
                xlabel('number of sessions')
                ylabel('RMS')
                axis square
                drawnow
%}

                % check signal differences between original & inverted
                orgpP.A = BPA.Ep.A;
                y    = spm_int_J(orgpP,M,U);            % integrate with observer

                j = ceil(ySt/2):floor(T/2);
                spm_figure('Getwin','Figure 3'); clf
                subplot(2,1,1)
                plot(t(j),orgY(j,:))
                subplot(2,1,2)
                plot(t(j),y(j,:))
                xlabel('Time (seconds)')
                ylabel('Amplitude')

                % MAE of signals
                mae = mean(mean(abs(y-orgY)))
                MAEs(k,i) = mae;

                % show original & simulated signal correlation        
                figure; Corr(k,1) = plotTwoSignalsCorrelation(y, orgY);
            catch exception
                msgText = getReport(exception)
            end

            % save inversion result
            save(fname, 'pPs', 'Uus', 'BPAs', 'RMSs', 'MAEs', 'invTms', 'CSD', 'Corr');
        end

    end % for k=1:maxK
end
