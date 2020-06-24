%==========================================================================
% Simulate timeseries
% Before running this script, please add path of spm12 directory
%
function reproductCheckSpm12DCMLargefMRI
    rng('default')

    % DEM Structure: create random inputs
    % -------------------------------------------------------------------------
    N  = 12;                              % number of runs
    T  = 200;                             % number of observations (scans)
    TR = 2;                               % repetition time or timing
    n  = 8;                               % number of regions or nodes
    t  = (1:T)*TR;                        % observation times

    % priors
    % -------------------------------------------------------------------------
    options.maxnodes   = n;               % effective number of nodes

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

    % integrate states
    % -------------------------------------------------------------------------
    U.dt = TR;
    M.f  = 'spm_fx_fmri';
    M.x  = sparse(n,5);
    M.g  = 'spm_gx_fmri';

    maxK = 20; % number of trial

    % load spectrum DCM performance check result
    fname = ['results/DCM_demo-rand' num2str(n) '-' num2str(N) 'x' num2str(maxK) '.mat'];
    load(fname);

    % finding minimum error row
    maeMin = min(min(MAEs));
    idx = find(MAEs==maeMin);
    k = mod(idx-1,50) + 1;
    minI = floor((idx-1)/50 + 1);

    ySt = 30; % BOLD signal starting point (for DCM inversion)
    T = 199 + ySt;

    %pP.A = rand(n,n)/4 - 0.125;
    pP = pPs{k}; % A,B,C,D params of minimum error row
    orgpP = pP;

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

    % init for inverted result
    CSD   = {};

    for i = 1:N
        try
            % integrate states
            % -----------------------------------------------------------------
            U.u  = Uus{k,i};                     % endogenous fluctuations
            y    = spm_int_J(pP,M,U);            % integrate with observer
            orgY = y;

            % response
            % -----------------------------------------------------------------
            DCM.Y.y  = y(ySt:end,:);

            % nonlinear system identification (Variational Laplace)
            % =================================================================
            CSD{end + 1} = spm_dcm_fmri_csd(DCM);
            BPA          = spm_dcm_average(CSD,'simulation',1);

            % summary
            % -----------------------------------------------------------------
            spm_figure('Getwin','Figure 2'); clf

            subplot(2,1,1); hold off
            spm_plot_ci(BPA.Ep.A(:),BPA.Cp(1:n*n,1:n*n)), hold on
            bar(pP.A(:),1/4), hold off
            title('True and MAP connections','FontSize',16)
            axis square

            % check signal differences between original & inverted
            orgpP.A = BPA.Ep.A;
            y = spm_int_J(orgpP,M,U);            % integrate with observer

            j = ceil(ySt/2):floor(T/2);
            spm_figure('Getwin','Figure 3'); clf
            subplot(2,1,1)
            plot(t(j),orgY(j,:))
            xlabel('Time (seconds)')
            ylabel('Amplitude')
            subplot(2,1,2)
            plot(t(j),y(j,:))
            xlabel('Time (seconds)')
            ylabel('Amplitude')

            % MAE of signals
            mae = mean(mean(abs(y-orgY)))

            % plot original A
            figure; plotDcmEC(pP.A);
            % plot predicted A
            figure; plotDcmEC(BPA.Ep.A);

        catch exception
            msgText = getReport(exception)
        end
    end
end
