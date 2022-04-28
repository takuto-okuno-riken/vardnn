% Before using this function, download PartiallyConditionedGrangerCausality codes from
% https://github.com/danielemarinazzo/PartiallyConditionedGrangerCausality
% and add a path "PartiallyConditionedGrangerCausality-master" and sub folders. 

function performanceCheckNodePatternTVB2
    node_nums = [11,22,33,44,55,66];
    num_scan = 55;
    if num_scan == 47 % deco's 66 node. weight add. VARDNN-GC show highest AUC, others so so.
        Gths = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
    elseif num_scan == 48 % deco's 66 node. original. FC so so. GC and VARDNN-GC show low AUC
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 2048 % deco's 66 node. original TR=2 BOLD / FC ok. others bad.
        Gths = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
    elseif num_scan == 3048 % deco's 66 node. original TR=0.1 BOLD / FC and VARDNN-GC so so.
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
    pcpcAUC = zeros(1,N);
    lsopcAUC = zeros(1,N);
    plspcAUC = zeros(1,N);
    wcsAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    dlmAUC = zeros(1,N);
    pcdlAUC = zeros(1,N);
    pcdlwAUC = zeros(1,N);
    dlgAUC = zeros(1,N);
    linueAUC = zeros(1,N);
    pcsAUC = zeros(1,N);
    cpcAUC = zeros(1,N);
    fgesAUC = zeros(1,N);
    fcaAUC = zeros(1,N);
    tsfcAUC = zeros(1,N);
    tsfcaAUC = zeros(1,N);
    mvardiAUC = zeros(1,N);
    pvardiAUC = zeros(1,N);
    mpcvardiAUC = zeros(1,N);
    mpcvargcAUC = zeros(1,N);
    ppcvardiAUC = zeros(1,N);
    ppcvargcAUC = zeros(1,N);
    mplsdiAUC = zeros(1,N);
    mplsgcAUC = zeros(1,N);
    pplsdiAUC = zeros(1,N);
    pplsgcAUC = zeros(1,N);
    mlsodiAUC = zeros(1,N);
    mlsogcAUC = zeros(1,N);
    pcgcAUC = zeros(1,N);
    dls1AUC = zeros(1,N);
    dls3AUC = zeros(1,N);
    msvmdiAUC = zeros(1,N);
    msvmgcAUC = zeros(1,N);
    mgpdiAUC = zeros(1,N);
    mgpgcAUC = zeros(1,N);
    mgpediAUC = zeros(1,N);
    nvdiAUC = zeros(1,N);
    nvmiAUC = zeros(1,N);
    trdiAUC = zeros(1,N);
    trmiAUC = zeros(1,N);
    rfdiAUC = zeros(1,N);
    rfmiAUC = zeros(1,N);
    svlpcAUC = zeros(1,N);
    svgpcAUC = zeros(1,N);
    svrpcAUC = zeros(1,N);
    gppcAUC = zeros(1,N);
    trpcAUC = zeros(1,N);
    rfpcAUC = zeros(1,N);
    ccmAUC = zeros(1,N);
    ccmpgAUC = zeros(1,N);
    ccmmgAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    pcpcROC = cell(N,2);
    lsopcROC = cell(N,2);
    plspcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);
    dlmROC = cell(N,2);
    pcdlROC = cell(N,2);
    pcdlwROC = cell(N,2);
    dlgROC = cell(N,2);
    linueROC = cell(N,2);
    pcsROC = cell(N,2);
    cpcROC = cell(N,2);
    fgesROC = cell(N,2);
    fcaROC = cell(N,2);
    tsfcROC = cell(N,2);
    tsfcaROC = cell(N,2);
    mvardiROC = cell(N,2);
    pvardiROC = cell(N,2);
    mpcvardiROC = cell(N,2);
    mpcvargcROC = cell(N,2);
    ppcvardiROC = cell(N,2);
    ppcvargcROC = cell(N,2);
    mplsdiROC = cell(N,2);
    mplsgcROC = cell(N,2);
    pplsdiROC = cell(N,2);
    pplsgcROC = cell(N,2);
    mlsodiROC = cell(N,2);
    mlsogcROC = cell(N,2);
    pcgcROC = cell(N,2);
    dls1ROC = cell(N,2);
    dls3ROC = cell(N,2);
    msvmdiROC = cell(N,2);
    msvmgcROC = cell(N,2);
    mgpdiROC = cell(N,2);
    mgpgcROC = cell(N,2);
    mgpediROC = cell(N,2);
    nvdiROC = cell(N,2);
    nvmiROC = cell(N,2);
    trdiROC = cell(N,2);
    trmiROC = cell(N,2);
    rfdiROC = cell(N,2);
    rfmiROC = cell(N,2);
    svlpcROC = cell(N,2);
    svgpcROC = cell(N,2);
    svrpcROC = cell(N,2);
    gppcROC = cell(N,2);
    trpcROC = cell(N,2);
    rfpcROC = cell(N,2);
    ccmROC = cell(N,2);
    ccmpgROC = cell(N,2);
    ccmmgROC = cell(N,2);
    fcRf = figure;
    pcRf = figure;
    pcpcRf = figure;
    lsopcRf = figure;
    plspcRf = figure;
    wcsRf = figure;
    gcRf = figure;
    pgcRf = figure;
    dlRf = figure;
    dlwRf = figure;
    dlmRf = figure;
    pcdlRf = figure;
    pcdlwRf = figure;
    dlgRf = figure;
    linueRf = figure;
    pcsRf = figure;
    cpcRf = figure;
    fgesRf = figure;
    fcaRf = figure;
    tsfcRf = figure;
    tsfcaRf = figure;
    mvardiRf = figure;
    pvardiRf = figure;
    mpcvardiRf = figure;
    mpcvargcRf = figure;
    ppcvardiRf = figure;
    ppcvargcRf = figure;
    mplsdiRf = figure;
    mplsgcRf = figure;
    pplsdiRf = figure;
    pplsgcRf = figure;
    mlsodiRf = figure;
    mlsogcRf = figure;
    pcgcRf = figure;
    dls1Rf = figure;
    dls3Rf = figure;
    msvmdiRf = figure;
    msvmgcRf = figure;
    mgpdiRf = figure;
    mgpgcRf = figure;
    mgpediRf = figure;
    nvdiRf = figure;
    nvmiRf = figure;
    trdiRf = figure;
    trmiRf = figure;
    rfdiRf = figure;
    rfmiRf = figure;
    svlpcRf = figure;
    svgpcRf = figure;
    svrpcRf = figure;
    gppcRf = figure;
    trpcRf = figure;
    rfpcRf = figure;
    ccmRf = figure;
    ccmpgRf = figure;
    ccmmgRf = figure;

    origf = figure;
    origSigf = figure;

    lag = 3;
    resfname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-result.mat'];
    if exist(resfname, 'file')
        load(resfname);
    end

    for k=1:N
        tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
        load(tvbFile);
        density = length(find(weights>Gth)) / (node_num * (node_num-1));
        nodeNum = size(si,1);
        sigLen = size(si,2);

        % show original connection
        figure(origf); plotDirectedFC(weights, 'Ground Truth', 1);
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
        if isempty(plspcROC{k,1})
            PC2 = calcPLSPartialCorrelation(si); % calc PLS PC
            Z = PC - PC2; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-PLSPC=' num2str(pcdiff)]);
    %        figure; clims = [-1 1]; imagesc(Z,clims); title('PC - PLSPC');
            figure(plspcRf); hold on; [plspcROC{k,1}, plspcROC{k,2}, plspcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title(['ROC curve of PC (pat=' num2str(i) ')']);
        end
        % show result of PCA-PC
        if isempty(pcpcROC{k,1})
            PC2 = calcPcPartialCorrelation(si); % calc PCA PC
            Z = PC - PC2; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-PCAPC=' num2str(pcdiff)]);
            figure(pcpcRf); hold on; [pcpcROC{k,1}, pcpcROC{k,2}, pcpcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title(['ROC curve of PCA-PC (pat=' num2str(i) ')']);
        end
        % show result of LassoPC
        if isempty(lsopcROC{k,1})
            [lambda, alpha, errMat] = estimateLassoParamsForPC(si, [], [], [], 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
            PC2 = calcLassoPartialCorrelation(si, [], [], [], lambda, alpha); % calc Lasso PC
            Z = PC - PC2; pcdiff=nanmean(abs(Z),'all'); disp(['mae of PC-LassoPC=' num2str(pcdiff)]);
            figure(lsopcRf); hold on; [lsopcROC{k,1}, lsopcROC{k,2}, lsopcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title(['ROC curve of Lasso PC (pat=' num2str(i) ')']);
        end
        % show result of SVR-PC
        if isempty(svlpcROC{k,1})
            PC2 = calcSvPartialCorrelation(si,[], [], [], 'linear');
            figure(svlpcRf); hold on; [svlpcROC{k,1}, svlpcROC{k,2}, svlpcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title('SVl-PC');
        end
        if isempty(svgpcROC{k,1})
            PC2 = calcSvPartialCorrelation(si,[], [], [], 'gaussian');
            figure(svgpcRf); hold on; [svgpcROC{k,1}, svgpcROC{k,2}, svgpcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title('SVg-PC');
        end
        if isempty(svrpcROC{k,1})
            PC2 = calcSvPartialCorrelation(si,[], [], [], 'rbf');
            figure(svrpcRf); hold on; [svrpcROC{k,1}, svrpcROC{k,2}, svrpcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title('SVr-PC');
        end
        % show result of GP-PC
        if isempty(gppcROC{k,1})
            PC2 = calcGpPartialCorrelation(si);
            figure(gppcRf); hold on; [gppcROC{k,1}, gppcROC{k,2}, gppcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title('GP-PC');
        end
        % show result of Tree-PC
        if isempty(trpcROC{k,1})
            PC2 = calcTreePartialCorrelation(si);
            figure(trpcRf); hold on; [trpcROC{k,1}, trpcROC{k,2}, trpcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title('Tree-PC');
        end
        % show result of RF-PC
        if isempty(rfpcROC{k,1})
            PC2 = calcRfPartialCorrelation(si);
            figure(rfpcRf); hold on; [rfpcROC{k,1}, rfpcROC{k,2}, rfpcAUC(k)] = plotROCcurve(PC2, weights, 100, 1, Gth); hold off;
            title('RF-PC');
        end

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
        if isempty(gcROC{k,1})
            gcI = calcMultivariateGCI_(si,[],[],[],lag);
            figure(gcRf); hold on; [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mvGC (pat=' num2str(i) ')']);
        end
        % show result of granger causality index (pwGC)
        if isempty(pgcROC{k,1})
            pgcI = calcPairwiseGCI(si,[],[],[],lag);
            figure(pgcRf); hold on; [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = plotROCcurve(pgcI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of pwGC (pat=' num2str(i) ')']);
        end
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

        % calcurate and show VARDNN-GC
        dlGC = [];
        exControl = eye(nodeNum, nodeNum);
        netFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
            if exist('inSignal','var'), exSignal=inSignal; end % for compatibility
        else
            % train VARDNN    
            Y = si;
            exSignal = uu;
            % layer parameters
            weightFunc = @estimateInitWeightRoughHe;
            weightParam = [10];
            bias = 0.5;
            netDLCM = initMvarDnnNetwork(Y, exSignal, [], exControl); % weightFunc, weightParam, bias);
            % training VARDNN network
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
            [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);

            % recoverty training
            %[netDLCM, time] = recoveryTrainMvarDnnNetwork(Y, exSignal, [], exControl, netDLCM, options);
            save(netFile, 'netDLCM', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
        end
        if isempty(dlGC)
            % show VARDNN-GC
            dlGC = calcMvarDnnGCI(Y, exSignal, [], exControl, netDLCM);
            save(netFile, 'netDLCM', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2', 'dlGC');
        end
        
        % calc ROC curve
        if isempty(dlROC{k,1})
            figure(dlRf); hold on; [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of VARDNN-GC (pat=' num2str(i) ')']);
        end
        % show result of VARDNN weight causality index (VARDNN-WCI) as VARDNN-DI
        if isempty(dlwROC{k,1})
            fg = figure; dlw = plotMvarDnnDI(netDLCM, [], exControl); close(fg);
            figure(dlwRf); hold on; [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = plotROCcurve(dlw, weights, 100, 1, Gth); hold off;
            title(['ROC curve of VARDNN-WCI (pat=' num2str(i) ')']);
        end
        % show result of VARDNN mean impact value (VARDNN-MIV)
        if isempty(dlmROC{k,1})
            fg = figure; [~,dlm] = calcMvarDnnMIV(Y, exSignal, [], exControl, netDLCM); close(fg);
            figure(dlmRf); hold on; [dlmROC{k,1}, dlmROC{k,2}, dlmAUC(k)] = plotROCcurve(dlm, weights, 100, 1, Gth); hold off;
            title(['ROC curve of VARDNN-MIV (pat=' num2str(i) ')']);
        end
%%}
        % calcurate and show PC-VARDNN
        netFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) 'mpcvar.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            % train VARDNN    
            Y = si;
            exSignal = uu;
            % layer parameters
            netMPC = initMpcvarDnnNetwork(Y, exSignal, [], exControl); % weightFunc, weightParam, bias);
            % training PC-VARDNN network
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
            netMPC = trainMpcvarDnnNetwork(Y, exSignal, [], exControl, netMPC, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(netMPC);
            disp(['end training : rsme=' num2str(rsme)]);

            % recoverty training
            save(netFile, 'netMPC', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
        end
        % show result of PC-VARDNN-GC
        if isempty(pcdlROC{k,1})
            fg = figure; dlGC = plotMpcvarDnnGCI(Y, exSignal, [], exControl, netMPC, 0); close(fg);
            figure(pcdlRf); hold on; [pcdlROC{k,1}, pcdlROC{k,2}, pcdlAUC(k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of PC-VARDNN-GC (pat=' num2str(i) ')']);
        end
        % show result of PC-VARDNN-DI
        if isempty(pcdlwROC{k,1})
            fg = figure; dlw = plotMpcvarDnnDI(netMPC, [], exControl); close(fg);
            figure(pcdlwRf); hold on; [pcdlwROC{k,1}, pcdlwROC{k,2}, pcdlwAUC(k)] = plotROCcurve(dlw, weights, 100, 1, Gth); hold off;
            title(['ROC curve of PC-VARDNN-DI (pat=' num2str(i) ')']);
        end
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
        if isempty(fcaROC{k,1})
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
        end
        % show result of multivaliate Vector Auto-Regression DI
        if isempty(mvardiROC{k,1})
            netMVAR = initMvarNetwork(si, [], [], [], lag);
            mvarDI = calcMvarDI(netMVAR, [], []);
            figure(mvardiRf); hold on; [mvardiROC{k,1}, mvardiROC{k,2}, mvardiAUC(k)] = plotROCcurve(mvarDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mVAR-DI (pat=' num2str(i) ')']);

            % show result of pairwise Vector Auto-Regression DI
            pvarDI = calcPvarDI(si, [], [], [], lag);
            figure(pvardiRf); hold on; [pvardiROC{k,1}, pvardiROC{k,2}, pvardiAUC(k)] = plotROCcurve(pvarDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of pVAR-DI (pat=' num2str(i) ')']);
        end
        % show result of multivaliate PC Vector Auto-Regression DI
        if isempty(mpcvardiROC{k,1})
            netMPCVAR = initMpcvarNetwork(si, [], [], [], lag);
            mpcvarDI = calcMpcvarDI(netMPCVAR, [], []);
            figure(mpcvardiRf); hold on; [mpcvardiROC{k,1}, mpcvardiROC{k,2}, mpcvardiAUC(k)] = plotROCcurve(mpcvarDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mPCVAR-DI (pat=' num2str(i) ')']);

            % show result of multivaliate PC Vector Auto-Regression GC
            mpcvarGC = calcMpcvarGCI(si, [], [], [], netMPCVAR);
            figure(mpcvargcRf); hold on; [mpcvargcROC{k,1}, mpcvargcROC{k,2}, mpcvargcAUC(k)] = plotROCcurve(mpcvarGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mPCVAR-GC (pat=' num2str(i) ')']);
        end
        % show result of pairwise PC Vector Auto-Regression DI
        if isempty(ppcvardiROC{k,1})
            netPPCVAR = initPpcvarNetwork(si, [], [], [], lag);
            ppcvarDI = calcPpcvarDI(netPPCVAR, [], []);
            figure(ppcvardiRf); hold on; [ppcvardiROC{k,1}, ppcvardiROC{k,2}, ppcvardiAUC(k)] = plotROCcurve(ppcvarDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of pPCVAR-DI (pat=' num2str(i) ')']);

            % show result of pairwise PC Vector Auto-Regression GC
            ppcvarGC = calcPpcvarGCI(si, [], [], [], netPPCVAR);
            figure(ppcvargcRf); hold on; [ppcvargcROC{k,1}, ppcvargcROC{k,2}, ppcvargcAUC(k)] = plotROCcurve(ppcvarGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of pPCVAR-GC (pat=' num2str(i) ')']);
        end
        % show result of multivaliate PLS Vector Auto-Regression DI
        if isempty(mplsdiROC{k,1})
            netMPLSVAR = initMplsvarNetwork(si, [], [], [], lag);
            mplsvarDI = calcMplsvarDI(netMPLSVAR, [], []);
            figure(mplsdiRf); hold on; [mplsdiROC{k,1}, mplsdiROC{k,2}, mplsdiAUC(k)] = plotROCcurve(mplsvarDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mPLSVAR-DI (pat=' num2str(i) ')']);

            % show result of multivaliate PLS Vector Auto-Regression GC
            mplsvarGC = calcMplsvarGCI(si, [], [], [], netMPLSVAR);
            figure(mplsgcRf); hold on; [mplsgcROC{k,1}, mplsgcROC{k,2}, mplsgcAUC(k)] = plotROCcurve(mplsvarGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mPLSVAR-GC (pat=' num2str(i) ')']);

            % mvGC and mPLSVARGC diff check
            Z = gcI - mplsvarGC; gcdiff=nanmean(abs(Z),'all'); disp(['mae of mvGC-mplsvarGC=' num2str(gcdiff) ' / ' num2str(max(max(gcI)))]);
            figure; clims = [-1 1]; imagesc(Z,clims); title('mvGC - mplsvarGC');
        end
        % show result of pairwise PLS Vector Auto-Regression DI
        if isempty(pplsdiROC{k,1})
            netPPLSVAR = initPplsvarNetwork(si, [], [], [], lag);
            pplsvarDI = calcPplsvarDI(netPPLSVAR, [], []);
            figure(pplsdiRf); hold on; [pplsdiROC{k,1}, pplsdiROC{k,2}, pplsdiAUC(k)] = plotROCcurve(pplsvarDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of pPLSVAR-DI (pat=' num2str(i) ')']);

            % show result of pairwise PLS Vector Auto-Regression GC
            pplsvarGC = calcPplsvarGCI(si, [], [], [], netPPLSVAR);
            figure(pplsgcRf); hold on; [pplsgcROC{k,1}, pplsgcROC{k,2}, pplsgcAUC(k)] = plotROCcurve(pplsvarGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of pPLSVAR-GC (pat=' num2str(i) ')']);
        end
        % ---------
        % (multivaliate Lasso Vector Auto-Regression DI) without exogenous signals
        if isempty(mlsodiROC{k,1})
            [lambda, elaAlpha, errMat] = estimateLassoParamsForMvar(si, [], [], [], lag, 0.5, 5, [0.01:0.02:0.99],[1:-0.1:0.1]);
            netMVAR = initMlassovarNetwork(si, [], [], [], lag, lambda, elaAlpha);
            [mlsoDI,~,mlso] = calcMlassovarDI(netMVAR, [], []);
            figure(mlsodiRf); hold on; [mlsodiROC{k,1}, mlsodiROC{k,2}, mlsodiAUC(k)] = plotROCcurve(mlsoDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mLSOVAR-DI (pat=' num2str(i) ')']);

            % (multivaliate Lasso Vector Auto-Regression GC) without exogenous signals
            mlsoGC = calcMlassovarGCI(si, [], [], [], netMVAR);
            figure(mlsogcRf); hold on; [mlsogcROC{k,1}, mlsogcROC{k,2}, mlsogcAUC(k)] = plotROCcurve(mlsoGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mLSOVAR-GC (pat=' num2str(i) ')']);
        end
        % ---------
        if isempty(pcgcROC{k,1})
            pcGC = calcPCGC(si, lag, 0.8);
            figure(pcgcRf); hold on; [pcgcROC{k,1}, pcgcROC{k,2}, pcgcAUC(k)] = plotROCcurve(pcGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of PC-GC (pat=' num2str(i) ')']);
        end
        % calcurate and show VARLSTM-GC
        netFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) 'ls1.mat'];
        if exist(netFile, 'file')
            load(netFile);
            if exist('inSignal','var'), exSignal=inSignal; end % for compatibility
        else
            % train VARLSTM
            Y = si;
            exSignal = uu;
            exControl = eye(nodeNum, nodeNum);
            % layer parameters
            net1 = initMvarLstmNetwork(Y, exSignal, [], exControl,1);
            net3 = initMvarLstmNetwork(Y, exSignal, [], exControl,3);
            % training PC-VARLSTM network
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
            net1 = trainMvarLstmNetwork(Y, exSignal, [], exControl, net1, options);
            net3 = trainMvarLstmNetwork(Y, exSignal, [], exControl, net3, options);
            save(netFile, 'net1', 'net3', 'Y', 'exSignal', 'exControl');
        end
        % calc ROC curve
        if isempty(dls1ROC{k,1})
            dlGC = calcMvarLstmGCI(Y, exSignal, [], exControl, net1);
            figure(dls1Rf); hold on; [dls1ROC{k,1}, dls1ROC{k,2}, dls1AUC(k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of VARLSTM(1)-GC (pat=' num2str(i) ')']);
        end
        % calc ROC curve
        if isempty(dls3ROC{k,1})
            dlGC = calcMvarLstmGCI(Y, exSignal, [], exControl, net3);
            figure(dls3Rf); hold on; [dls3ROC{k,1}, dls3ROC{k,2}, dls3AUC(k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of VARLSTM(1)-GC (pat=' num2str(i) ')']);
        end
        % ---------
        % (multivaliate SVM Vector Auto-Regression DI) without exogenous signals
        if isempty(msvmdiROC{k,1})
            netMVAR = initMsvmvarNetwork(si, [], [], [], lag);
            [msvmDI] = calcMsvmvarDI(netMVAR, [], []);
            figure(msvmdiRf); hold on; [msvmdiROC{k,1}, msvmdiROC{k,2}, msvmdiAUC(k)] = plotROCcurve(msvmDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mSVMVAR-DI (pat=' num2str(i) ')']);

            % (multivaliate SVM Vector Auto-Regression GC) without exogenous signals
            msvmGC = calcMsvmvarGCI(si, [], [], [], netMVAR);
            figure(msvmgcRf); hold on; [msvmgcROC{k,1}, msvmgcROC{k,2}, msvmgcAUC(k)] = plotROCcurve(msvmGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mSVMVAR-GC (pat=' num2str(i) ')']);
        end
        % (multivaliate GP Vector Auto-Regression DI) without exogenous signals
        if isempty(mgpdiROC{k,1})
            netMVAR = initMgpvarNetwork(si, [], [], [], lag);
            [mgpDI] = calcMgpvarDI(netMVAR, [], []);
            figure(mgpdiRf); hold on; [mgpdiROC{k,1}, mgpdiROC{k,2}, mgpdiAUC(k)] = plotROCcurve(mgpDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mGPVAR-DI (pat=' num2str(i) ')']);

            % (multivaliate GP Vector Auto-Regression GC) without exogenous signals
            mgpGC = calcMgpvarGCI(si, [], [], [], netMVAR);
            figure(mgpgcRf); hold on; [mgpgcROC{k,1}, mgpgcROC{k,2}, mgpgcAUC(k)] = plotROCcurve(mgpGC, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mGPVAR-GC (pat=' num2str(i) ')']);

            netMVAR = initMgpvarNetwork(si, [], [], [], lag, 'ardsquaredexponential', 'constant');
            [mgpeDI] = calcMgpvarDI(netMVAR, [], []);
            figure(mgpediRf); hold on; [mgpediROC{k,1}, mgpediROC{k,2}, mgpediAUC(k)] = plotROCcurve(mgpeDI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of mGPeVAR-DI (pat=' num2str(i) ')']);
        end
        % extra tests (nVARNN DI)
        netFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) 'nvnn.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            exControl = ones(1,nodeNum);
            % layer parameters
            nvNN = initNvarnnNetwork(si, exSignal, [], exControl,1);
            % training NVARNN network
            maxEpochs = 2000;
            miniBatchSize = ceil(size(si,2) / 3);

            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',0.01, ...
                'Verbose',true);

            disp('start training');
            nvNN = trainNvarnnNetwork(si, exSignal, [], exControl, nvNN, options);
            save(netFile, 'si','exSignal','exControl','nvNN');
        end
        if isempty(nvdiROC{k,1})
            fg = figure; varDI = plotNvarnnDI(nvNN, [], exControl); close(fg);
            figure(nvdiRf); hold on; [nvdiROC{k,1}, nvdiROC{k,2}, nvdiAUC(k)] = plotROCcurve(varDI, weights, 100, 1, Gth); hold off;
            title('nVARNN-DI');
            % extra tests (nVARNN MIV)
            [~,varMAIV] = calcNvarnnMIV(si, exSignal,  [], exControl, nvNN);
            figure(nvmiRf); hold on; [nvmiROC{k,1}, nvmiROC{k,2}, nvmiAUC(k)] = plotROCcurve(varMAIV, weights, 100, 1, Gth); hold off;
            title('nVARNN-MAIV');
        end
        % extra tests (mTreeVAR DI)
        if isempty(trdiROC{k,1})
            netMVAR = initMtreevarNetwork(si, [], [], [], lag);
            fg = figure; varDI = plotMtreevarDI(netMVAR, [], []); close(fg);
            figure(trdiRf); hold on; [trdiROC{k,1}, trdiROC{k,2}, trdiAUC(k)] = plotROCcurve(varDI, weights, 100, 1, Gth); hold off;
            title('mTreeVAR-DI');
            % extra tests (mTreeVAR MIV)
            [~,trvarMAIV] = calcMtreevarMIV(si, [], [], [], netMVAR);
            figure(trmiRf); hold on; [trmiROC{k,1}, trmiROC{k,2}, trmiAUC(k)] = plotROCcurve(trvarMAIV, weights, 100, 1, Gth); hold off;
            title('mTreeVAR-MAIV');
        end
        % extra tests (mRFVAR DI)
        if isempty(rfdiROC{k,1})
            netMVAR = initMrfvarNetwork(si, [], [], [], lag);
            fg = figure; varDI = plotMrfvarDI(netMVAR, [], []); close(fg);
            figure(rfdiRf); hold on; [rfdiROC{k,1}, rfdiROC{k,2}, rfdiAUC(k)] = plotROCcurve(varDI, weights, 100, 1, Gth); hold off;
            title('mRFVAR-DI');
            % extra tests (mRFVAR MIV)
            [~,varMAIV] = calcMrfvarMIV(si, [], [], [], netMVAR);
            figure(rfmiRf); hold on; [rfmiROC{k,1}, rfmiROC{k,2}, rfmiAUC(k)] = plotROCcurve(varMAIV, weights, 100, 1, Gth); hold off;
            title('mRFVAR-MAIV');
        end
        % extra tests (CCM)
        if isempty(ccmROC{k,1})
            fg = figure; CCM = plotConvCrossMap(si, [], [], [], lag); close(fg);
            figure(ccmRf); hold on; [ccmROC{k,1}, ccmROC{k,2}, ccmAUC(k)] = plotROCcurve(CCM, weights, 100, 1, Gth); hold off;
            title('pwCCM');
        end
        % extra tests (CCM pGC & mGC)
        if isempty(ccmpgROC{k,1})
            % CCM-pGC
            fg = figure; CCM = plotConvCrossMapPGC(si, [], [], [], lag); close(fg);
            figure(ccmpgRf); hold on; [ccmpgROC{k,1}, ccmpgROC{k,2}, ccmpgAUC(k)] = plotROCcurve(CCM, weights, 100, 1, Gth); hold off;
            title('CCM-pGC');
            % CCM-mGC
            fg = figure; CCM = plotConvCrossMapMGC(si, [], [], [], lag); close(fg);
            figure(ccmmgRf); hold on; [ccmmgROC{k,1}, ccmmgROC{k,2}, ccmmgAUC(k)] = plotROCcurve(CCM, weights, 100, 1, Gth); hold off;
            title('CCM-mGC');
        end
    end
    % show result AUC
    disp(['FC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(fcAUC))]);
    disp(['PC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pcAUC))]);
    disp(['PCA-PC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pcpcAUC))]);
    disp(['WCS AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(wcsAUC))]);
    disp(['mvGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(gcAUC))]);
    disp(['pwGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pgcAUC))]);
    disp(['VARDNN-GC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlAUC))]);
    disp(['VARDNN-WCI AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlwAUC))]);
    disp(['VARDNN-MIV AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlmAUC))]);
    disp(['LINUE-TE AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(linueAUC))]);
    disp(['DirectLiNGAM AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlgAUC))]);
    disp(['PC-sm AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pcsAUC))]);

    % save result
    save(resfname, 'fcAUC','pcAUC','pcpcAUC','lsopcAUC','plspcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','dlmAUC','pcdlAUC','pcdlwAUC','dlgAUC','linueAUC','pcsAUC','cpcAUC','fgesAUC','fcaAUC','tsfcAUC','tsfcaAUC', ...
        'mvardiAUC','mpcvardiAUC','mpcvargcAUC','pvardiAUC','ppcvardiAUC','ppcvargcAUC','mplsdiAUC','mplsgcAUC','pplsdiAUC','pplsgcAUC','mlsodiAUC','mlsogcAUC','pcgcAUC','dls1AUC','dls3AUC','msvmdiAUC','msvmgcAUC','mgpdiAUC','mgpgcAUC','mgpediAUC', ...
        'nvdiAUC','nvmiAUC','trdiAUC','trmiAUC','rfdiAUC','rfmiAUC', ...
        'svlpcAUC','svgpcAUC','svrpcAUC','gppcAUC','trpcAUC','rfpcAUC', ...
        'ccmAUC','ccmpgAUC','ccmmgAUC', ...
        'fcROC','pcROC','pcpcROC','lsopcROC','plspcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlmROC','pcdlROC','pcdlwROC','dlgROC','linueROC','pcsROC','cpcROC','fgesROC','fcaROC','tsfcROC','tsfcaROC', ...
        'mvardiROC','mpcvardiROC','mpcvargcROC','pvardiROC','ppcvardiROC','ppcvargcROC','mplsdiROC','mplsgcROC','pplsdiROC','pplsgcROC','mlsodiROC','mlsogcROC','pcgcROC','dls1ROC','dls3ROC','msvmdiROC','msvmgcROC','mgpdiROC','mgpgcROC','mgpediROC', ...
        'nvdiROC','nvmiROC','trdiROC','trmiROC','rfdiROC','rfmiROC', ...
        'svlpcROC','svgpcROC','svrpcROC','gppcROC','trpcROC','rfpcROC', ...
        'ccmROC','ccmpgROC','ccmmgROC');

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
    plotAverageROCcurve(mvardiROC, N, '-', [0.3,0.3,0.3],0.5);
    plotAverageROCcurve(pvardiROC, N, '--', [0.3,0.3,0.3],0.5);
    plotAverageROCcurve(mpcvardiROC, N, '-', [0.3,0.6,0.6],1.0);
    plotAverageROCcurve(mpcvargcROC, N, '--', [0.3,0.6,0.6],0.8);
    plotAverageROCcurve(ppcvardiROC, N, '-', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(ppcvargcROC, N, '--', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(mplsdiROC, N, '-', [0.7,0.9,0.9],1.0);
    plotAverageROCcurve(mplsgcROC, N, '--', [0.7,0.9,0.9],0.8);
    plotAverageROCcurve(pplsdiROC, N, '-', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(pplsgcROC, N, '--', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(mlsodiROC, N, '-.', [0.7,0.9,0.9],1.0);
    plotAverageROCcurve(mlsogcROC, N, ':', [0.7,0.9,0.9],0.8);
    plotAverageROCcurve(pcgcROC, N, '-.', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(ccmROC, N, '--', [0.9,0.2,0.2],0.8);
    plotAverageROCcurve(ccmpgROC, N, '-.', [0.9,0.2,0.2],0.8);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(i)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
