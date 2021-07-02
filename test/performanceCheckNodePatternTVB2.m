function performanceCheckNodePatternTVB2
    node_nums = [11,22,33,44,55,66];
    num_scan = 55;
    if num_scan == 47 % deco's 66 node. weight add. DLCM-GC show highest AUC, others so so.
        Gths = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
    elseif num_scan == 48 % deco's 66 node. original. FC so so. GC and DLCM-GC show low AUC
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 2048 % deco's 66 node. original TR=2 BOLD / FC ok. others bad.
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 3048 % deco's 66 node. original TR=0.1 BOLD / FC and DLCM-GC so so.
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 2948 % deco's 66 node. original TR=0.1 BOLD / all bad.
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 3148 % deco's 66 node. original TR=0.2 BOLD / all bad.
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 50  % TVB's 66 node. 
        Gths = [0.02, 0.01, 0.01, 0.01, 0.01, 0.01];
    elseif num_scan == 51  % oh's mouse 98 node. original.
        node_nums = [16,32,48,64,80,98];
        Gths = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05];
    elseif num_scan == 52  % oh's mouse 98 node. density around 0.15. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
    elseif num_scan == 53  % oh's mouse 98 node. density around 0.23. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
    elseif num_scan == 54  % oh's mouse 98 node. density 0.15. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [1, 1, 1, 1, 1, 1];
    elseif num_scan == 55  % oh's mouse 98 node. density 0.15. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [1, 1, 1, 1, 1, 1];
    end
    % test sparse and full density
    hz = 64;
    N = 8;

    for i=1:length(node_nums)
        checkingPattern(node_nums(i), num_scan, hz, Gths(i), N, i);
    end
end

function checkingPattern(node_num, num_scan, hz, Gth, N, i)
    % init
    fcAUC = zeros(1,N);
    pcAUC = zeros(1,N);
    wcsAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    dlgAUC = zeros(1,N);
    linueAUC = zeros(1,N);
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
    mplsvarecAUC = zeros(1,N);
    mplsvargcAUC = zeros(1,N);
    pplsvarecAUC = zeros(1,N);
    pplsvargcAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);
    dlgROC = cell(N,2);
    linueROC = cell(N,2);
    pcsROC = cell(N,2);
    cpcROC = cell(N,2);
    fgesROC = cell(N,2);
    fcaROC = cell(N,2);
    tsfcROC = cell(N,2);
    tsfcaROC = cell(N,2);
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
    fcRf = figure;
    pcRf = figure;
    wcsRf = figure;
    gcRf = figure;
    pgcRf = figure;
    dlRf = figure;
    dlwRf = figure;
    dlgRf = figure;
    linueRf = figure;
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
    mplsvarecRf = figure;
    mplsvargcRf = figure;
    pplsvarecRf = figure;
    pplsvargcRf = figure;

    origf = figure;
    origSigf = figure;

    lag = 3;

    for k=1:N
        tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
        load(tvbFile);
        density = length(find(weights>Gth)) / (node_num * (node_num-1));
        nodeNum = size(si,1);
        sigLen = size(si,2);

        % show original connection
        figure(origf); plotEC(weights, 'Ground Truth', 1);
        figure(origSigf); plot(t, si);

        % show result of FC
        FC = calcFunctionalConnectivity(si);
        figure(fcRf); hold on; [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = plotROCcurve(FC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of FC (pat=' num2str(i) ')']);

        % show result of PC
        PC = calcPartialCorrelation(si);
        figure(pcRf); hold on; [pcROC{k,1}, pcROC{k,2}, pcAUC(k)] = plotROCcurve(PC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of PC (pat=' num2str(i) ')']);
        % PC and PLSPC diff check
        PLSPC = calcPLSPartialCorrelation(si); % calc PLS PC
        Z = PC - PLSPC; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-PLSPC=' num2str(pcdiff)]);
        figure; clims = [-1 1]; imagesc(Z,clims); title('PC - PLSPC');

        % show result of WCS
        wcsFile = ['results/wcs-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(wcsFile, 'file')
            load(wcsFile);
        else
            WCS = calcWaveletCoherence(si);
            save(wcsFile, 'WCS');
        end
        figure(wcsRf); hold on; [wcsROC{k,1}, wcsROC{k,2}, wcsAUC(k)] = plotROCcurve(WCS, weights, 100, 1, Gth); hold off;
        title(['ROC curve of WCS (pat=' num2str(i) ')']);

        % show result of granger causality index (mvGC)
        gcI = calcMultivariateGCI_(si,[],[],[],lag);
        figure(gcRf); hold on; [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mvGC (pat=' num2str(i) ')']);

        % show result of granger causality index (pwGC)
        pgcI = calcPairwiseGCI(si,[],[],[],lag);
        figure(pgcRf); hold on; [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = plotROCcurve(pgcI, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pwGC (pat=' num2str(i) ')']);

        % show result of DirectLiNGAM
        dlgFile = ['results/dling-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(dlgFile, 'file')
            load(dlgFile);
        else
            Aest = calcDirectLiNGAM(si);
            save(dlgFile, 'Aest');
        end
        figure(dlgRf); hold on; [dlgROC{k,1}, dlgROC{k,2}, dlgAUC(k)] = plotROCcurve(Aest, weights, 100, 1, Gth); hold off;
        title(['ROC curve of DirectLiNGAM (pat=' num2str(i) ')']);
%%{
        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(si);
        [uu, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            
        % calcurate and show DLCM-GC
        dlGC = [];
        exControl = eye(nodeNum, nodeNum);
        dlcmFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
            if exist('inSignal','var'), exSignal=inSignal; end % for compatibility
        else
            % train DLCM    
            Y = si;
            exSignal = uu;
            % layer parameters
            weightFunc = @estimateInitWeightRoughHe;
            weightParam = [10];
            bias = 0.5;
            netDLCM = initMvarDnnNetwork(Y, exSignal, [], exControl); % weightFunc, weightParam, bias);
            % training DLCM network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 3);
            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.05, ...
                'Verbose',false);
        %            'Plots','training-progress');

            disp('start training');
            netDLCM = trainMvarDnnNetwork(Y, exSignal, [], exControl, netDLCM, options);
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);

            % recoverty training
            %[netDLCM, time] = recoveryTrainDlcmNetwork(Y, exSignal, [], exControl, netDLCM, options);
            save(dlcmFile, 'netDLCM', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
        end
        if isempty(dlGC)
            % show DLCM-GC
            dlGC = calcMvarDnnGCI(Y, exSignal, [], exControl, netDLCM);
            save(dlcmFile, 'netDLCM', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2', 'dlGC');
        end
        
        % calc ROC curve
        figure(dlRf); hold on; [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of DLCM-GC (pat=' num2str(i) ')']);

        % show result of DLCM weight causality index (DLCM-wci) as DLCM-EC
        fg = figure; dlwGC = plotMvarDnnEC(netDLCM, [], exControl); close(fg);
        figure(dlwRf); hold on; [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = plotROCcurve(dlwGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of DLCM-WCI (pat=' num2str(i) ')']);
%%}
%%{
        % linue TE result
        linueFile = ['results/tvb-pat-linue/linue_MultivAnalysis_tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '-' num2str(lag) '.mat'];
        load(linueFile);
        A = outputToStore.reshapedMtx.';

        % show ROC curve of TE(LIN UE)
        figure(linueRf); hold on; [linueROC{k,1}, linueROC{k,2}, linueAUC(k)] = plotROCcurve(A, weights, 100, 1, Gth); hold off;        
        title(['ROC curve of LINUE-TE (pat=' num2str(i) ')']);
%%}
        % show result of TETRAD PC-stable-max
        csvFile = ['results/tetrad/pcs-tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.csv'];
        A = readmatrix(csvFile);
        figure(pcsRf); hold on; [pcsROC{k,1}, pcsROC{k,2}, pcsAUC(k)] = plotROCcurve(A, weights, 100, 1, Gth); hold off;
        title(['ROC curve of PC-stable-max (pat=' num2str(i) ')']);

        % show result of TETRAD CPC
        csvFile = ['results/tetrad/cpc-tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.csv'];
        A = readmatrix(csvFile);
        figure(cpcRf); hold on; [cpcROC{k,1}, cpcROC{k,2}, cpcAUC(k)] = plotROCcurve(A, weights, 100, 1, Gth); hold off;
        title(['ROC curve of CPC (pat=' num2str(i) ')']);

        % show result of TETRAD FGES
        csvFile = ['results/tetrad/fges-tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.csv'];
        A = readmatrix(csvFile);
        figure(fgesRf); hold on; [fgesROC{k,1}, fgesROC{k,2}, fgesAUC(k)] = plotROCcurve(A, weights, 100, 1, Gth); hold off;
        title(['ROC curve of CPC (pat=' num2str(i) ')']);
        
        % -----------------------------------------------------------------
        % extra tests
        % show result of FC abs
        FCa = calcFunctionalConnectivityAbs(si);
        figure(fcaRf); hold on; [fcaROC{k,1}, fcaROC{k,2}, fcaAUC(k)] = plotROCcurve(FCa, weights, 100, 1, Gth); hold off;
        title(['ROC curve of FCa (pat=' num2str(i) ')']);

        % show result of time shifted FC
        tsFC = calcTimeShiftedCorrelation(si,[],[],[],lag);
        figure(tsfcRf); hold on; [tsfcROC{k,1}, tsfcROC{k,2}, tsfcAUC(k)] = plotROCcurve(tsFC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of tsFC (pat=' num2str(i) ')']);

        % show result of time shifted FC abs
        tsFCa = calcTimeShiftedCorrelationAbs(si,[],[],[],lag);
        figure(tsfcaRf); hold on; [tsfcaROC{k,1}, tsfcaROC{k,2}, tsfcaAUC(k)] = plotROCcurve(tsFCa, weights, 100, 1, Gth); hold off;
        title(['ROC curve of tsFCa (pat=' num2str(i) ')']);

        % show result of multivaliate Vector Auto-Regression EC
        netMVAR = initMvarNetwork(si, [], [], [], lag);
        mvarEC = calcMvarEC(netMVAR, [], []);
        figure(mvarecRf); hold on; [mvarecROC{k,1}, mvarecROC{k,2}, mvarecAUC(k)] = plotROCcurve(mvarEC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mVAR-EC (pat=' num2str(i) ')']);

        % show result of pairwise Vector Auto-Regression EC
        pvarEC = calcPvarEC(si, [], [], [], lag);
        figure(pvarecRf); hold on; [pvarecROC{k,1}, pvarecROC{k,2}, pvarecAUC(k)] = plotROCcurve(pvarEC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pVAR-EC (pat=' num2str(i) ')']);

        % show result of multivaliate PC Vector Auto-Regression EC
        netMPCVAR = initMpcvarNetwork(si, [], [], [], lag);
        mpcvarEC = calcMpcvarEC(netMPCVAR, [], []);
        figure(mpcvarecRf); hold on; [mpcvarecROC{k,1}, mpcvarecROC{k,2}, mpcvarecAUC(k)] = plotROCcurve(mpcvarEC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mPCVAR-EC (pat=' num2str(i) ')']);

        % show result of multivaliate PC Vector Auto-Regression GC
        mpcvarGC = calcMpcvarGCI(si, [], [], [], netMPCVAR);
        figure(mpcvargcRf); hold on; [mpcvargcROC{k,1}, mpcvargcROC{k,2}, mpcvargcAUC(k)] = plotROCcurve(mpcvarGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mPCVAR-GC (pat=' num2str(i) ')']);

        % show result of pairwise PC Vector Auto-Regression EC
        netPPCVAR = initPpcvarNetwork(si, [], [], [], lag);
        ppcvarEC = calcPpcvarEC(netPPCVAR, [], []);
        figure(ppcvarecRf); hold on; [ppcvarecROC{k,1}, ppcvarecROC{k,2}, ppcvarecAUC(k)] = plotROCcurve(ppcvarEC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pPCVAR-EC (pat=' num2str(i) ')']);

        % show result of pairwise PC Vector Auto-Regression GC
        ppcvarGC = calcPpcvarGCI(si, [], [], [], netPPCVAR);
        figure(ppcvargcRf); hold on; [ppcvargcROC{k,1}, ppcvargcROC{k,2}, ppcvargcAUC(k)] = plotROCcurve(ppcvarGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pPCVAR-GC (pat=' num2str(i) ')']);

        % show result of multivaliate PLS Vector Auto-Regression EC
        netMPLSVAR = initMplsvarNetwork(si, [], [], [], lag);
        mplsvarEC = calcMplsvarEC(netMPLSVAR, [], []);
        figure(mplsvarecRf); hold on; [mplsvarecROC{k,1}, mplsvarecROC{k,2}, mplsvarecAUC(k)] = plotROCcurve(mplsvarEC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mPLSVAR-EC (pat=' num2str(i) ')']);

        % show result of multivaliate PLS Vector Auto-Regression GC
        mplsvarGC = calcMplsvarGCI(si, [], [], [], netMPLSVAR);
        figure(mplsvargcRf); hold on; [mplsvargcROC{k,1}, mplsvargcROC{k,2}, mplsvargcAUC(k)] = plotROCcurve(mplsvarGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mPLSVAR-GC (pat=' num2str(i) ')']);

        % mvGC and mPLSVARGC diff check
        Z = gcI - mplsvarGC; gcdiff=nanmean(abs(Z),'all'); disp(['mae of mvGC-mplsvarGC=' num2str(gcdiff) ' / ' num2str(max(max(gcI)))]);
        figure; clims = [-1 1]; imagesc(Z,clims); title('mvGC - mplsvarGC');

        % show result of pairwise PLS Vector Auto-Regression EC
        netPPLSVAR = initPplsvarNetwork(si, [], [], [], lag);
        pplsvarEC = calcPplsvarEC(netPPLSVAR, [], []);
        figure(pplsvarecRf); hold on; [pplsvarecROC{k,1}, pplsvarecROC{k,2}, pplsvarecAUC(k)] = plotROCcurve(pplsvarEC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pPLSVAR-EC (pat=' num2str(i) ')']);

        % show result of pairwise PLS Vector Auto-Regression GC
        pplsvarGC = calcPplsvarGCI(si, [], [], [], netPPLSVAR);
        figure(pplsvargcRf); hold on; [pplsvargcROC{k,1}, pplsvargcROC{k,2}, pplsvargcAUC(k)] = plotROCcurve(pplsvarGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pPLSVAR-GC (pat=' num2str(i) ')']);
    end
    % show result AUC
    disp(['FC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(fcAUC))]);
    disp(['PC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pcAUC))]);
    disp(['WCS AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(wcsAUC))]);
    disp(['mvGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(gcAUC))]);
    disp(['pwGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pgcAUC))]);
    disp(['DLCM-GC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlAUC))]);
    disp(['DLCM-WCI AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlwAUC))]);
    disp(['LINUE-TE AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(linueAUC))]);
    disp(['DirectLiNGAM AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlgAUC))]);
    disp(['PC-sm AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pcsAUC))]);

    % save result
    fname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-result.mat'];
    save(fname, 'fcAUC','pcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','dlgAUC','linueAUC','pcsAUC','cpcAUC','fgesAUC','fcaAUC','tsfcAUC','tsfcaAUC','mvarecAUC','mpcvarecAUC', ...
        'fcROC','pcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlgROC','linueROC','pcsROC','cpcROC','fgesROC','fcaROC','tsfcROC','tsfcaROC','mvarecROC','mpcvarecROC');

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
%    plotErrorROCcurve(dcmROC, N, [0.2,0.2,0.8]);
%    plotErrorROCcurve(rnnROC, N, [0.8,0.8,0.2]);
    plotErrorROCcurve(linueROC, N, [0.2,0.6,0.8]);
%    plotErrorROCcurve(nnnueROC, N, [0.8,0.2,0.8]);
    plotErrorROCcurve(dlgROC, N, [0.6,0.6,0.3]);
    plotErrorROCcurve(pcsROC, N, [0.5,0.5,0.5]);
    plotErrorROCcurve(cpcROC, N, [0.5,0.5,0.5]);
    plotErrorROCcurve(fgesROC, N, [0.5,0.5,0.5]);
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '--', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '--', [0.0,0.5,0.0],0.5);
%    plotAverageROCcurve(dlROC, N, '-', [0.2,0.2,0.2],1.2);
%    plotAverageROCcurve(dlwROC, N, '--', [0.2,0.2,0.2],0.8); % TODO:
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.2);
    plotAverageROCcurve(dlROC, N, '--', [0.2,0.2,0.2],0.8); % TODO:
%    plotAverageROCcurve(dcmROC, N, '-', [0.2,0.2,0.8],0.5);
%    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5);
    plotAverageROCcurve(linueROC, N, '--', [0.2,0.5,0.8],0.5);
%    plotAverageROCcurve(nnnueROC, N, '--', [0.7,0.2,0.7],0.5);
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
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(i)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
