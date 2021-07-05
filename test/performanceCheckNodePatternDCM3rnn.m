% Before using this function, download SPM12 codes from
% https://www.fil.ion.ucl.ac.uk/spm/software/download/
% and add a path "spm12" and sub folders, then remove "spm12/external" folder and sub folders.

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
    plspcAUC = zeros(1,N);
    mvarecAUC = zeros(1,N);
    pvarecAUC = zeros(1,N);
    mpcvarecAUC = zeros(1,N);
    mpcvargcAUC = zeros(1,N);
    ppcvarecAUC = zeros(1,N);
    ppcvargcAUC = zeros(1,N);
    mplsvarecAUC = zeros(1,N);
    mplsvargcAUC = zeros(1,N);
    pplsvarecAUC = zeros(1,N);
    pplsvargcAUC = zeros(1,N);
    plsoecAUC = zeros(1,N);
    mlsoecAUC = zeros(1,N);
    plsogcAUC = zeros(1,N);
    mlsogcAUC = zeros(1,N);
    rnnROC = cell(N,2);
    linueROC = cell(N,2);
    nnnueROC = cell(N,2);
    pcsROC = cell(N,2);
    fgesROC = cell(N,2);
    cpcROC = cell(N,2);
    fcaROC = cell(N,2);
    tsfcROC = cell(N,2);
    tsfcaROC = cell(N,2);
    plspcROC = cell(N,2);
    mvarecROC = cell(N,2);
    pvarecROC = cell(N,2);
    mpcvarecROC = cell(N,2);
    mpcvargcROC = cell(N,2);
    ppcvarecROC = cell(N,2);
    ppcvargcROC = cell(N,2);
    mplsvarecROC = cell(N,2);
    mplsvargcROC = cell(N,2);
    pplsvarecROC = cell(N,2);
    pplsvargcROC = cell(N,2);
    plsoecROC = cell(N,2);
    mlsoecROC = cell(N,2);
    plsogcROC = cell(N,2);
    mlsogcROC = cell(N,2);
    rnnRf = figure;
    linueRf = figure;
    nnnueRf = figure;
    pcsRf = figure;
    cpcRf = figure;
    fgesRf = figure;
    fcaRf = figure;
    tsfcRf = figure;
    tsfcaRf = figure;
    plspcRf = figure;
    mvarecRf = figure;
    pvarecRf = figure;
    mpcvarecRf = figure;
    mpcvargcRf = figure;
    ppcvarecRf = figure;
    ppcvargcRf = figure;
    mplsvarecRf = figure;
    mplsvargcRf = figure;
    pplsvarecRf = figure;
    pplsvargcRf = figure;
    plsoecRf = figure;
    mlsoecRf = figure;
    plsogcRf = figure;
    mlsogcRf = figure;
    
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
        % extra tests (PLS PC)
        fg = figure; PLSPC = plotFunctionalConnectivityAbs(y2.'); close(fg);
        figure(plspcRf); hold on; [plspcROC{k,1}, plspcROC{k,2}, plspcAUC(k)] = plotROCcurve(PLSPC, pP.A, 100, 1, Gth); hold off;
        title('PLS PC');
        % extra tests (multivaliate Vector Auto-Regression EC)
        netMVAR = initMvarNetwork(y2.', [], [], [], 3);
        fg = figure; mvarEC = plotMvarEC(netMVAR, [], []); close(fg);
        figure(mvarecRf); hold on; [mvarecROC{k,1}, mvarecROC{k,2}, mvarecAUC(k)] = plotROCcurve(mvarEC, pP.A, 100, 1, Gth); hold off;
        title('mVAR-EC');
        % extra tests (pairwised Vector Auto-Regression EC)
        fg = figure; pvarEC = plotPvarEC(y2.', [], [], [], 3); close(fg);
        figure(pvarecRf); hold on; [pvarecROC{k,1}, pvarecROC{k,2}, pvarecAUC(k)] = plotROCcurve(pvarEC, pP.A, 100, 1, Gth); hold off;
        title('pVAR-EC');
        % extra tests (multivaliate Principal Component Vector Auto-Regression EC)
        netMPCVAR = initMpcvarNetwork(y2.', [], [], [], 3);
        fg = figure; mpcvarEC = plotMpcvarEC(netMPCVAR, [], []); close(fg);
        figure(mpcvarecRf); hold on; [mpcvarecROC{k,1}, mpcvarecROC{k,2}, mpcvarecAUC(k)] = plotROCcurve(mpcvarEC, pP.A, 100, 1, Gth); hold off;
        title('mPCVAR-EC');
        % extra tests (multivaliate Principal Component Vector Auto-Regression GC)
        fg = figure; mpcvarGC = plotMpcvarGCI(y2.', [], [], [], netMPCVAR); close(fg);
        figure(mpcvargcRf); hold on; [mpcvargcROC{k,1}, mpcvargcROC{k,2}, mpcvargcAUC(k)] = plotROCcurve(mpcvarGC, pP.A, 100, 1, Gth); hold off;
        title('mPCVAR-GC');
        % extra tests (pairwise Principal Component Vector Auto-Regression EC)
        netPPCVAR = initPpcvarNetwork(y2.', [], [], [], 3);
        fg = figure; ppcvarEC = plotPpcvarEC(netPPCVAR, [], []); close(fg);
        figure(ppcvarecRf); hold on; [ppcvarecROC{k,1}, ppcvarecROC{k,2}, ppcvarecAUC(k)] = plotROCcurve(ppcvarEC, pP.A, 100, 1, Gth); hold off;
        title('pPCVAR-EC');
        % extra tests (pairwise Principal Component Vector Auto-Regression GC)
        fg = figure; ppcvarGC = plotPpcvarGCI(y2.', [], [], [], netPPCVAR); close(fg);
        figure(ppcvargcRf); hold on; [ppcvargcROC{k,1}, ppcvargcROC{k,2}, ppcvargcAUC(k)] = plotROCcurve(ppcvarGC, pP.A, 100, 1, Gth); hold off;
        title('pPCVAR-GC');
        % extra tests (multivaliate PLS Vector Auto-Regression EC)
        netMPLSVAR = initMplsvarNetwork(y2.', [], [], [], 3);
        fg = figure; mplsvarEC = plotMplsvarEC(netMPLSVAR, [], []); close(fg);
        figure(mplsvarecRf); hold on; [mplsvarecROC{k,1}, mplsvarecROC{k,2}, mplsvarecAUC(k)] = plotROCcurve(mplsvarEC, pP.A, 100, 1, Gth); hold off;
        title('mPLSVAR-EC');
        % extra tests (multivaliate PLS Vector Auto-Regression GC)
        fg = figure; mplsvarGC = plotMplsvarGCI(y2.', [], [], [], netMPLSVAR); close(fg);
        figure(mplsvargcRf); hold on; [mplsvargcROC{k,1}, mplsvargcROC{k,2}, mplsvargcAUC(k)] = plotROCcurve(mplsvarGC, pP.A, 100, 1, Gth); hold off;
        title('mPLSVAR-GC');
        % extra tests (pairwise PLS Vector Auto-Regression EC)
        netPPLSVAR = initPplsvarNetwork(y2.', [], [], [], 3);
        fg = figure; pplsvarEC = plotPplsvarEC(netPPLSVAR, [], []); close(fg);
        figure(pplsvarecRf); hold on; [pplsvarecROC{k,1}, pplsvarecROC{k,2}, pplsvarecAUC(k)] = plotROCcurve(pplsvarEC, pP.A, 100, 1, Gth); hold off;
        title('pPLSVAR-EC');
        % extra tests (pairwise PLS Vector Auto-Regression GC)
        fg = figure; pplsvarGC = plotPplsvarGCI(y2.', [], [], [], netPPLSVAR); close(fg);
        figure(pplsvargcRf); hold on; [pplsvargcROC{k,1}, pplsvargcROC{k,2}, pplsvargcAUC(k)] = plotROCcurve(pplsvarGC, pP.A, 100, 1, Gth); hold off;
        title('pPLSVAR-GC');
        % show result of pLassoVAR EC
        [lambda, elaAlpha, errMat] = estimateLassoParamsForMvar(y2.', [], [], [], 3, 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
        fg = figure; EC = plotPlassovarDI(y2.', [], [], [], 3, lambda, elaAlpha); close(fg);
        figure(plsoecRf); hold on; [plsoecROC{k,1}, plsoecROC{k,2}, plsoecAUC(k)] = plotROCcurve(EC, pP.A); hold off;
        title('pLassoVAR-EC');
        % show result of mLassoVAR EC
        netLsoVAR = initMlassovarNetwork(y2.', [], [], [], 3, lambda, elaAlpha);
        fg = figure; EC = plotMlassovarDI(netLsoVAR, [], []); close(fg);
        figure(mlsoecRf); hold on; [mlsoecROC{k,1}, mlsoecROC{k,2}, mlsoecAUC(k)] = plotROCcurve(EC, pP.A); hold off;
        title('mLassoVAR-EC');
        % show result of pLassoVAR GC
        fg = figure; GC = plotPlassovarGCI(y2.', [], [], [], 3, 0.01); close(fg);
        figure(plsogcRf); hold on; [plsogcROC{k,1}, plsogcROC{k,2}, plsogcAUC(k)] = plotROCcurve(GC, pP.A); hold off;
        title('pLassoVAR-GC');
        % show result of mLassoVAR GC
        fg = figure; GC = plotMlassovarGCI(y2.', [], [], [], netLsoVAR); close(fg);
        figure(mlsogcRf); hold on; [mlsogcROC{k,1}, mlsogcROC{k,2}, mlsogcAUC(k)] = plotROCcurve(GC, pP.A); hold off;
        title('mLassoVAR-GC');
    end
    % save result
    save(fname, 'fcAUC','pcAUC','pcpcAUC','plspcAUC','lsopcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','dlgAUC','dcmAUC', ...
        'rnnAUC','linueAUC','nnnueAUC','pcsAUC','cpcAUC','fgesAUC','fcaAUC','tsfcAUC','tsfcaAUC', ...
        'mvarecAUC','pvarecAUC','mpcvarecAUC','mpcvargcAUC','ppcvarecAUC','ppcvargcAUC',...
        'mplsvarecAUC','mplsvargcAUC','pplsvarecAUC','pplsvargcAUC',...
        'plsoecAUC','mlsoecAUC','plsogcAUC','mlsogcAUC',...
        'fcROC','pcROC','pcpcROC','plspcROC','lsopcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlgROC','dcmROC', ...
        'rnnROC','linueROC','nnnueROC','pcsROC','cpcROC','fgesROC','fcaROC','tsfcROC','tsfcaROC', ...
        'mvarecROC','pvarecROC','mpcvarecROC','mpcvargcROC','ppcvarecROC','ppcvargcROC', ...
        'mplsvarecROC','mplsvargcROC','pplsvarecROC','pplsvargcROC', ...
        'plsoecROC','mlsoecROC','plsogcROC','mlsogcROC');
    
    % show average ROC curve of DCM
    figure; 
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(wcsROC, N, [0.9,0.5,0]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(pgcROC, N, [0.0,0.5,0.0]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % TODO:
    plotErrorROCcurve(dcmROC, N, [0.2,0.2,0.8]);
    plotErrorROCcurve(rnnROC, N, [0.8,0.8,0.2]);
    plotErrorROCcurve(linueROC, N, [0.2,0.6,0.8]);
    plotErrorROCcurve(nnnueROC, N, [0.8,0.2,0.8]);
    plotErrorROCcurve(dlgROC, N, [0.6,0.6,0.3]);
    plotErrorROCcurve(pcsROC, N, [0.5,0.5,0.5]);
    plotErrorROCcurve(cpcROC, N, [0.5,0.5,0.5]);
    plotErrorROCcurve(fgesROC, N, [0.5,0.5,0.5]);
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '-', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pcpcROC, N, '--', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(lsopcROC, N, '-.', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(plspcROC, N, ':', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '--', [0.0,0.5,0.0],0.5);
%    plotAverageROCcurve(dlROC, N, '-', [0.2,0.2,0.2],1.2);
%    plotAverageROCcurve(dlwROC, N, '--', [0.2,0.2,0.2],0.8); % TODO:
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
    plotAverageROCcurve(mplsvarecROC, N, '-', [0.7,0.9,0.9],1.0);
    plotAverageROCcurve(mplsvargcROC, N, '--', [0.7,0.9,0.9],0.8);
    plotAverageROCcurve(pplsvarecROC, N, '-', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(pplsvargcROC, N, '--', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(plsoecROC, N, '-', [0.9,0.6,0.9],0.5);
    plotAverageROCcurve(plsogcROC, N, '--', [0.9,0.6,0.9],0.5);
    plotAverageROCcurve(mlsoecROC, N, '-', [0.9,0.7,0.9],1.0);
    plotAverageROCcurve(mlsogcROC, N, '--', [0.9,0.7,0.9],0.8);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(idx)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

function [x, y] = plotAverageROCcurve(roc, N, line, col, width)
    x = zeros(N,size(roc{1,1},2));
    y = zeros(N,size(roc{1,1},2));
    for k=1:N
        x(k,:) = roc{k,1};
        y(k,:) = roc{k,2};
    end
    mx = mean(x,1);
    my = mean(y,1);
    % plot line
    plot(mx, my, line,'Color',col,'LineWidth',width);
end

function [x, y] = plotErrorROCcurve(roc, N, col)
    x = zeros(N,size(roc{1,1},2));
    y = zeros(N,size(roc{1,1},2));
    for k=1:N
        x(k,:) = roc{k,1};
        y(k,:) = roc{k,2};
    end
    mx = mean(x,1);
    my = mean(y,1);
    % error area
    errx = std(x,1,1) / sqrt(size(x,1));
    erry = std(y,1,1) / sqrt(size(y,1));
    xt = mx -errx;
    yt = my +erry;
    xb = mx +errx;
    yb = my -erry;
    Y = [yb fliplr(yt)];
    X = [xb fliplr(xt)];
    patch(X,Y,col,'EdgeColor','none','FaceAlpha',.04);
end
