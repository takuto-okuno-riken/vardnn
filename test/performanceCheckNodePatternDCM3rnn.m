% Before using this function, download SPM12 codes from
% https://www.fil.ion.ucl.ac.uk/spm/software/download/
% and add a path "spm12" and sub folders, then remove "spm12/external" folder and sub folders.

% Before using this function, download PartiallyConditionedGrangerCausality codes from
% https://github.com/danielemarinazzo/PartiallyConditionedGrangerCausality
% and add a path "PartiallyConditionedGrangerCausality-master" and sub folders. 

% this script should run after performanceCheckNodePatternDCM3, performanceCheckNodePatternDCM3d and RNN-GC result
function performanceCheckNodePatternDCM3rnn
    % -------------------------------------------------------------------------
    N  = 8;
    T  = 300;                             % number of observations (scans)
    n  = 8;                               % number of regions or nodes

    prefix = 'net-pat3-';                 % original weight file prefix (result of *NodePatternDCM3d.m)
    Gth = 0;                              % 0 for pat3. 0.2 for pat4.
%    prefix = 'net-pat4-';                  % original weight file prefix (result of *NodePatternDCM3d.m)
%    Gth = 0.2;                             % 0 for pat3. 0.2 for pat4.

    %% pattern 1 -------------------------------------------------
%%{
    disp('network density 0.05');
    checkingPattern(N,T,n,prefix,Gth,1);
%%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('network density 0.11');
    checkingPattern(N,T,n,prefix,Gth,2);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('network density 0.16');
    checkingPattern(N,T,n,prefix,Gth,6);
%%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('network density 0.27');
    checkingPattern(N,T,n,prefix,Gth,7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('network density 0.36');
    checkingPattern(N,T,n,prefix,Gth,8);
%%}
end

%% 
function checkingPattern(N,T,n,prefix,Gth,idx)
    fname = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    load(fname);

    % init
    rnnAUC = zeros(1,N);
    linueAUC = zeros(1,N);
    nnnueAUC = zeros(1,N);
    pcsAUC = zeros(1,N);
    cpcAUC = zeros(1,N);
    fgesAUC = zeros(1,N);
    fcaAUC = zeros(1,N);
    tsfcAUC = zeros(1,N);
    tsfcaAUC = zeros(1,N);
    mvarecAUC = zeros(1,N);
    pvarecAUC = zeros(1,N);
    mpcvarecAUC = zeros(1,N);
    mpcvargcAUC = zeros(1,N);
    ppcvarecAUC = zeros(1,N);
    ppcvargcAUC = zeros(1,N);
    mplsecAUC = zeros(1,N);
    mplsgcAUC = zeros(1,N);
    pplsecAUC = zeros(1,N);
    pplsgcAUC = zeros(1,N);
    plsoecAUC = zeros(1,N);
    mlsoecAUC = zeros(1,N);
    plsogcAUC = zeros(1,N);
    mlsogcAUC = zeros(1,N);
    pcgcAUC = zeros(1,N);
    dls1AUC = zeros(1,N);
    dls3AUC = zeros(1,N);
    rnnROC = cell(N,2);
    linueROC = cell(N,2);
    nnnueROC = cell(N,2);
    pcsROC = cell(N,2);
    fgesROC = cell(N,2);
    cpcROC = cell(N,2);
    fcaROC = cell(N,2);
    tsfcROC = cell(N,2);
    tsfcaROC = cell(N,2);
    mvarecROC = cell(N,2);
    pvarecROC = cell(N,2);
    mpcvarecROC = cell(N,2);
    mpcvargcROC = cell(N,2);
    ppcvarecROC = cell(N,2);
    ppcvargcROC = cell(N,2);
    mplsecROC = cell(N,2);
    mplsgcROC = cell(N,2);
    pplsecROC = cell(N,2);
    pplsgcROC = cell(N,2);
    plsoecROC = cell(N,2);
    mlsoecROC = cell(N,2);
    plsogcROC = cell(N,2);
    mlsogcROC = cell(N,2);
    pcgcROC = cell(N,2);
    dls1ROC = cell(N,2);
    dls3ROC = cell(N,2);
    rnnRf = figure;
    linueRf = figure;
    nnnueRf = figure;
    pcsRf = figure;
    cpcRf = figure;
    fgesRf = figure;
    fcaRf = figure;
    tsfcRf = figure;
    tsfcaRf = figure;
    mvarecRf = figure;
    pvarecRf = figure;
    mpcvarecRf = figure;
    mpcvargcRf = figure;
    ppcvarecRf = figure;
    ppcvargcRf = figure;
    mplsecRf = figure;
    mplsgcRf = figure;
    pplsecRf = figure;
    pplsgcRf = figure;
    plsoecRf = figure;
    mlsoecRf = figure;
    plsogcRf = figure;
    mlsogcRf = figure;
    pcgcRf = figure;
    dls1Rf = figure;
    dls3Rf = figure;
    
    origf = figure;
    rnnTrial = 8;

    % reading RNN-GC, TE(LIN UE), TE(BIN NUE), TETRAD algorithms results
    for k=1:N
        netFile = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(netFile);

        % show original connection
        figure(origf); plotDcmEC(pP.A);

        % -----------------------------------------------------------------
        % read RNN-GC result
        A = zeros(n,n);
        for j=1:rnnTrial
            rnnFile = ['results/rnngc/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' sprintf('%02d',k) '_' num2str(j) '.txt'];
            Aj = readmatrix(rnnFile);
            A = A + Aj.';
            %figure; plotDcmEC(Aj.');
        end
        A = A / rnnTrial;
        %figure; plotDcmEC(A);

        % show ROC curve of RNN-GC
        figure(rnnRf); hold on; [rnnROC{k,1}, rnnROC{k,2}, rnnAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;
        
        % -----------------------------------------------------------------
        % read Transfer Entropy (LIN UE) result
        linueFile = ['results/linue/linue_MultivAnalysis_' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(linueFile);
        A = outputToStore.reshapedMtx.';

        % show ROC curve of TE(LIN UE)
        figure(linueRf); hold on; [linueROC{k,1}, linueROC{k,2}, linueAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;        

        % -----------------------------------------------------------------
        % read Transfer Entropy (NearestNeighber NUE) result
        nnnueFile = ['results/nnnue/nnnue_MultivAnalysis_' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(nnnueFile);
        A = outputToStore.reshapedMtx.';

        % show ROC curve of TE(NearestNeighber NUE)
        figure(nnnueRf); hold on; [nnnueROC{k,1}, nnnueROC{k,2}, nnnueAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;        
        
        % -----------------------------------------------------------------
        % read TETRAD PC-stable-max result
        csvFile = ['results/tetrad/pcs-' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.csv'];
        A = readmatrix(csvFile);

        % show ROC curve of PC-stable-max
        figure(pcsRf); hold on; [pcsROC{k,1}, pcsROC{k,2}, pcsAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;

        % -----------------------------------------------------------------
        % read TETRAD CPC result
        csvFile = ['results/tetrad/cpc-' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.csv'];
        A = readmatrix(csvFile);

        % show ROC curve of CPC result
        figure(cpcRf); hold on; [cpcROC{k,1}, cpcROC{k,2}, cpcAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;

        % -----------------------------------------------------------------
        % read TETRAD FGES result
        csvFile = ['results/tetrad/fges-' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.csv'];
        A = readmatrix(csvFile);

        % show ROC curve of FGES result
        figure(fgesRf); hold on; [fgesROC{k,1}, fgesROC{k,2}, fgesAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;

        % -----------------------------------------------------------------
        % extra tests (FC abs)
        if isempty(fcaROC{k,1})
            fg = figure; FCa = plotFunctionalConnectivityAbs(y2.'); close(fg);
            figure(fcaRf); hold on; [fcaROC{k,1}, fcaROC{k,2}, fcaAUC(k)] = plotROCcurve(FCa, pP.A, 100, 1, Gth); hold off;
            title('FCa');
            % show result of time shifted FC
            fg = figure; tsFC = plotTimeShiftedCorrelation(y2.', [], [], [], 3); close(fg);
            figure(tsfcRf); hold on; [tsfcROC{k,1}, tsfcROC{k,2}, tsfcAUC(k)] = plotROCcurve(tsFC, pP.A, 100, 1, Gth); hold off;
            title('tsFC');
            % show result of time shifted FC (abs)
            fg = figure; tsFCa = plotTimeShiftedCorrelationAbs(y2.', [], [], [], 3); close(fg);
            figure(tsfcaRf); hold on; [tsfcaROC{k,1}, tsfcaROC{k,2}, tsfcaAUC(k)] = plotROCcurve(tsFCa, pP.A, 100, 1, Gth); hold off;
            title('tsFCa');
        end
        % extra tests (multivaliate Vector Auto-Regression DI)
        if isempty(mvarecROC{k,1})
            netMVAR = initMvarNetwork(y2.', [], [], [], 3);
            fg = figure; mvarDI = plotMvarDI(netMVAR, [], []); close(fg);
            figure(mvarecRf); hold on; [mvarecROC{k,1}, mvarecROC{k,2}, mvarecAUC(k)] = plotROCcurve(mvarDI, pP.A, 100, 1, Gth); hold off;
            title('mVAR-DI');
            % extra tests (pairwised Vector Auto-Regression DI)
            fg = figure; pvarDI = plotPvarDI(y2.', [], [], [], 3); close(fg);
            figure(pvarecRf); hold on; [pvarecROC{k,1}, pvarecROC{k,2}, pvarecAUC(k)] = plotROCcurve(pvarDI, pP.A, 100, 1, Gth); hold off;
            title('pVAR-DI');
        end
        % extra tests (multivaliate Principal Component Vector Auto-Regression DI)
        if isempty(mpcvarecROC{k,1})
            netMPCVAR = initMpcvarNetwork(y2.', [], [], [], 3);
            fg = figure; mpcvarDI = plotMpcvarDI(netMPCVAR, [], []); close(fg);
            figure(mpcvarecRf); hold on; [mpcvarecROC{k,1}, mpcvarecROC{k,2}, mpcvarecAUC(k)] = plotROCcurve(mpcvarDI, pP.A, 100, 1, Gth); hold off;
            title('mPCVAR-DI');
            % extra tests (multivaliate Principal Component Vector Auto-Regression GC)
            fg = figure; mpcvarGC = plotMpcvarGCI(y2.', [], [], [], netMPCVAR); close(fg);
            figure(mpcvargcRf); hold on; [mpcvargcROC{k,1}, mpcvargcROC{k,2}, mpcvargcAUC(k)] = plotROCcurve(mpcvarGC, pP.A, 100, 1, Gth); hold off;
            title('mPCVAR-GC');
        end
        % extra tests (pairwise Principal Component Vector Auto-Regression DI)
        if isempty(ppcvarecROC{k,1})
            netPPCVAR = initPpcvarNetwork(y2.', [], [], [], 3);
            fg = figure; ppcvarDI = plotPpcvarDI(netPPCVAR, [], []); close(fg);
            figure(ppcvarecRf); hold on; [ppcvarecROC{k,1}, ppcvarecROC{k,2}, ppcvarecAUC(k)] = plotROCcurve(ppcvarDI, pP.A, 100, 1, Gth); hold off;
            title('pPCVAR-DI');
            % extra tests (pairwise Principal Component Vector Auto-Regression GC)
            fg = figure; ppcvarGC = plotPpcvarGCI(y2.', [], [], [], netPPCVAR); close(fg);
            figure(ppcvargcRf); hold on; [ppcvargcROC{k,1}, ppcvargcROC{k,2}, ppcvargcAUC(k)] = plotROCcurve(ppcvarGC, pP.A, 100, 1, Gth); hold off;
            title('pPCVAR-GC');
        end
        % extra tests (multivaliate PLS Vector Auto-Regression DI)
        if isempty(mplsecROC{k,1})
            netMPLSVAR = initMplsvarNetwork(y2.', [], [], [], 3);
            fg = figure; mplsvarDI = plotMplsvarDI(netMPLSVAR, [], []); close(fg);
            figure(mplsecRf); hold on; [mplsecROC{k,1}, mplsecROC{k,2}, mplsecAUC(k)] = plotROCcurve(mplsvarDI, pP.A, 100, 1, Gth); hold off;
            title('mPLSVAR-DI');
            % extra tests (multivaliate PLS Vector Auto-Regression GC)
            fg = figure; mplsvarGC = plotMplsvarGCI(y2.', [], [], [], netMPLSVAR); close(fg);
            figure(mplsgcRf); hold on; [mplsgcROC{k,1}, mplsgcROC{k,2}, mplsgcAUC(k)] = plotROCcurve(mplsvarGC, pP.A, 100, 1, Gth); hold off;
            title('mPLSVAR-GC');
        end
        % extra tests (pairwise PLS Vector Auto-Regression DI)
        if isempty(pplsecROC{k,1})
            netPPLSVAR = initPplsvarNetwork(y2.', [], [], [], 3);
            fg = figure; pplsvarDI = plotPplsvarDI(netPPLSVAR, [], []); close(fg);
            figure(pplsecRf); hold on; [pplsecROC{k,1}, pplsecROC{k,2}, pplsecAUC(k)] = plotROCcurve(pplsvarDI, pP.A, 100, 1, Gth); hold off;
            title('pPLSVAR-DI');
            % extra tests (pairwise PLS Vector Auto-Regression GC)
            fg = figure; pplsvarGC = plotPplsvarGCI(y2.', [], [], [], netPPLSVAR); close(fg);
            figure(pplsgcRf); hold on; [pplsgcROC{k,1}, pplsgcROC{k,2}, pplsgcAUC(k)] = plotROCcurve(pplsvarGC, pP.A, 100, 1, Gth); hold off;
            title('pPLSVAR-GC');
        end
        % show result of pLassoVAR DI
        if isempty(plsoecROC{k,1})
            lsoFile = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) 'lso.mat'];
            if exist(lsoFile, 'file')
                load(lsoFile);
            else
                [lambda, elaAlpha, errMat] = estimateLassoParamsForMvar(y2.', [], [], [], 3, 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
                save(lsoFile, 'lambda','elaAlpha','errMat');
            end
            fg = figure; DI = plotPlassovarDI(y2.', [], [], [], 3, lambda, elaAlpha); close(fg);
            figure(plsoecRf); hold on; [plsoecROC{k,1}, plsoecROC{k,2}, plsoecAUC(k)] = plotROCcurve(DI, pP.A, 100, 1, Gth); hold off;
            title('pLassoVAR-DI');
            % show result of mLassoVAR DI
            netLsoVAR = initMlassovarNetwork(y2.', [], [], [], 3, lambda, elaAlpha);
            fg = figure; DI = plotMlassovarDI(netLsoVAR, [], []); close(fg);
            figure(mlsoecRf); hold on; [mlsoecROC{k,1}, mlsoecROC{k,2}, mlsoecAUC(k)] = plotROCcurve(DI, pP.A, 100, 1, Gth); hold off;
            title('mLassoVAR-DI');
            % show result of pLassoVAR GC
            fg = figure; GC = plotPlassovarGCI(y2.', [], [], [], 3, 0.01); close(fg);
            figure(plsogcRf); hold on; [plsogcROC{k,1}, plsogcROC{k,2}, plsogcAUC(k)] = plotROCcurve(GC, pP.A, 100, 1, Gth); hold off;
            title('pLassoVAR-GC');
            % show result of mLassoVAR GC
            fg = figure; GC = plotMlassovarGCI(y2.', [], [], [], netLsoVAR); close(fg);
            figure(mlsogcRf); hold on; [mlsogcROC{k,1}, mlsogcROC{k,2}, mlsogcAUC(k)] = plotROCcurve(GC, pP.A, 100, 1, Gth); hold off;
            title('mLassoVAR-GC');
        end
        % show result of PCGC
        if isempty(pcgcROC{k,1})
            fg = figure; GC = plotPCGC(y2.', 3, 0.8, 0); close(fg);
            figure(pcgcRf); hold on; [pcgcROC{k,1}, pcgcROC{k,2}, pcgcAUC(k)] = plotROCcurve(GC, pP.A, 100, 1, Gth); hold off;
            title('PC-GC');
        end
        % VARLSTM-GC
        netFile = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) 'ls1.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(y2.', 0);
            [exSignal, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(u2.', 0);
            exControl = eye(n,n);
            % layer parameters
            net1 = initMvarLstmNetwork(si, exSignal, [], exControl,1);
            net3 = initMvarLstmNetwork(si, exSignal, [], exControl,3);
            % training PC-VARDNN network
            maxEpochs = 1000;
            miniBatchSize = ceil(size(si,2) / 3);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'SequenceLength','longest', ...
                'Shuffle','never', ...
                'GradientThreshold',5,...
                'Verbose',false);

            disp('start training');
            net1 = trainMvarLstmNetwork(si, exSignal, [], exControl, net1, options);
            net3 = trainMvarLstmNetwork(si, exSignal, [], exControl, net3, options);
            save(netFile, 'si','exSignal','exControl','net1', 'net3');
        end
        % show result of VARLSTM-GC
        if isempty(dls1ROC{k,1})
            fg = figure; dlGC = plotMvarLstmGCI(si, exSignal, [], exControl, net1, 0); close(fg);
            figure(dls1Rf); hold on; [dls1ROC{k,1}, dls1ROC{k,2}, dls1AUC(k)] = plotROCcurve(dlGC, pP.A, 100, 1, Gth); hold off;
            title('VARLSTM(1)-GC');
        end
%        if isempty(dls3ROC{k,1})
            fg = figure; dlGC = plotMvarLstmGCI(si, exSignal, [], exControl, net3, 0); close(fg);
            figure(dls3Rf); hold on; [dls3ROC{k,1}, dls3ROC{k,2}, dls3AUC(k)] = plotROCcurve(dlGC, pP.A, 100, 1, Gth); hold off;
            title('VARLSTM(3)-GC');
%        end
    end
    % save result
    save(fname, 'fcAUC','pcAUC','pcpcAUC','plspcAUC','lsopcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','dlmAUC','dlgAUC','pcdlAUC','pcdlwAUC','dcmAUC', ...
        'rnnAUC','linueAUC','nnnueAUC','pcsAUC','cpcAUC','fgesAUC','fcaAUC','tsfcAUC','tsfcaAUC', ...
        'mvarecAUC','pvarecAUC','mpcvarecAUC','mpcvargcAUC','ppcvarecAUC','ppcvargcAUC',...
        'mplsecAUC','mplsgcAUC','pplsecAUC','pplsgcAUC',...
        'plsoecAUC','mlsoecAUC','plsogcAUC','mlsogcAUC','pcgcAUC','dls1AUC','dls3AUC',...
        'fcROC','pcROC','pcpcROC','plspcROC','lsopcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlmROC','dlgROC','pcdlROC','pcdlwROC','dcmROC', ...
        'rnnROC','linueROC','nnnueROC','pcsROC','cpcROC','fgesROC','fcaROC','tsfcROC','tsfcaROC', ...
        'mvarecROC','pvarecROC','mpcvarecROC','mpcvargcROC','ppcvarecROC','ppcvargcROC', ...
        'mplsecROC','mplsgcROC','pplsecROC','pplsgcROC', ...
        'plsoecROC','mlsoecROC','plsogcROC','mlsogcROC','pcgcROC','dls1ROC','dls3ROC');

    % show all average ROC curves
    figure; 
    hold on;
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '-', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pcpcROC, N, '--', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(lsopcROC, N, '-.', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(plspcROC, N, ':', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '--', [0.0,0.5,0.0],0.5);
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.2);
    plotAverageROCcurve(dlROC, N, '--', [0.2,0.2,0.2],0.8); % TODO:
    plotAverageROCcurve(dcmROC, N, '-', [0.2,0.2,0.8],0.5);
    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5);
    plotAverageROCcurve(linueROC, N, '--', [0.2,0.5,0.8],0.5);
    plotAverageROCcurve(nnnueROC, N, '--', [0.7,0.2,0.7],0.5);
    plotAverageROCcurve(dlgROC, N, '-.', [0.6,0.6,0.3],0.5);
    plotAverageROCcurve(pcsROC, N, '-', [0.5,0.5,0.5],0.5);
    plotAverageROCcurve(cpcROC, N, '--', [0.5,0.5,0.5],0.5);
    plotAverageROCcurve(fgesROC, N, '-.', [0.5,0.5,0.5],0.5);
%    plotAverageROCcurve(fcaROC, N, '-.', [0.8,0.2,0.2],0.5);
%    plotAverageROCcurve(tsfcROC, N, '-', [0.6,0.2,0.2],1.2);
%    plotAverageROCcurve(tsfcaROC, N, '-.', [0.6,0.2,0.2],1.2);
    plotAverageROCcurve(mvarecROC, N, '-', [0.3,0.3,0.3],0.5);
    plotAverageROCcurve(pvarecROC, N, '--', [0.3,0.3,0.3],0.5);
    plotAverageROCcurve(mpcvarecROC, N, '-', [0.3,0.6,0.6],1.0);
    plotAverageROCcurve(mpcvargcROC, N, '--', [0.3,0.6,0.6],0.8);
    plotAverageROCcurve(ppcvarecROC, N, '-', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(ppcvargcROC, N, '--', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(mplsecROC, N, '-', [0.7,0.9,0.9],1.0);
    plotAverageROCcurve(mplsgcROC, N, '--', [0.7,0.9,0.9],0.8);
    plotAverageROCcurve(pplsecROC, N, '-', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(pplsgcROC, N, '--', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(plsoecROC, N, '-', [0.9,0.6,0.9],0.5);
    plotAverageROCcurve(plsogcROC, N, '--', [0.9,0.6,0.9],0.5);
    plotAverageROCcurve(mlsoecROC, N, '-', [0.9,0.7,0.9],1.0);
    plotAverageROCcurve(mlsogcROC, N, '--', [0.9,0.7,0.9],0.8);
    plotAverageROCcurve(pcgcROC, N, '-.', [0.3,0.6,0.6],0.5);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(idx)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

