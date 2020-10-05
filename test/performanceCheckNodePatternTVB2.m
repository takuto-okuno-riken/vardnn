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
    dlgAUC = zeros(1,N);
    linueAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlgROC = cell(N,2);
    linueROC = cell(N,2);
    fcRf = figure;
    pcRf = figure;
    wcsRf = figure;
    gcRf = figure;
    pgcRf = figure;
    dlRf = figure;
    dlgRf = figure;
    linueRf = figure;
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
        figure(origf); plotDcmEC(weights);
        figure(origSigf); plot(t, si);

        % show original signal FC
        FC = calcFunctionalConnectivity(si);
        figure(fcRf); hold on; [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = plotROCcurve(FC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of FC (pat=' num2str(i) ')']);

        % show original signal PC
        PC = calcPartialCorrelation(si);
        figure(pcRf); hold on; [pcROC{k,1}, pcROC{k,2}, pcAUC(k)] = plotROCcurve(PC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of PC (pat=' num2str(i) ')']);

        % show original signal WCS
        wcsFile = ['results/wcs-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(wcsFile, 'file')
            load(wcsFile);
        else
            WCS = calcWaveletCoherence(si);
            save(wcsFile, 'WCS');
        end
        figure(wcsRf); hold on; [wcsROC{k,1}, wcsROC{k,2}, wcsAUC(k)] = plotROCcurve(WCS, weights, 100, 1, Gth); hold off;
        title(['ROC curve of WCS (pat=' num2str(i) ')']);

        % show original signal granger causality index (mvGC)
        gcI = calcMultivariateGCI(si,lag);
        figure(gcRf); hold on; [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mvGC (pat=' num2str(i) ')']);

        % show original signal granger causality index (pwGC)
        gcI = calcPairwiseGCI(si,lag);
        figure(pgcRf); hold on; [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pwGC (pat=' num2str(i) ')']);

        % show original signal DirectLiNGAM
        Aest = calcDirectLiNGAM(si);
        figure(dlgRf); hold on; [dlgROC{k,1}, dlgROC{k,2}, dlgAUC(k)] = plotROCcurve(Aest, weights, 100, 1, Gth); hold off;
        title(['ROC curve of DirectLiNGAM (pat=' num2str(i) ')']);
%%{
        [si, sig, m, maxsi, minsi] = convert2SigmoidSignal(si);
        [uu, sig2, m2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            
        % calcurate and show DLCM-GC
        inControl = eye(nodeNum, nodeNum);
        dlcmFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            % train DLCM    
            Y = si;
            inSignal = uu;
            % layer parameters
            weightFunc = @estimateInitWeightRoughHe;
            weightParam = [10];
            bias = 0.5;
            netDLCM = initDlcmNetwork(Y, inSignal, [], inControl); % weightFunc, weightParam, bias);
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
            netDLCM = trainDlcmNetwork(Y, inSignal, [], inControl, netDLCM, options);
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);

            % recoverty training
            %[netDLCM, time] = recoveryTrainDlcmNetwork(Y, inSignal, [], inControl, netDLCM, options);
            save(dlcmFile, 'netDLCM', 'Y', 'inSignal', 'si', 'sig', 'm', 'maxsi', 'minsi', 'sig2', 'm2', 'maxsi2', 'minsi2');
        end
        % show DLCM-GC
        dlGC = calcDlcmGCI(Y, inSignal, [], inControl, netDLCM);
        
        % calc ROC curve
        figure(dlRf); hold on; [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = plotROCcurve(dlGC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of DLCM-GC (pat=' num2str(i) ')']);
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
    end
    % show result AUC
    disp(['FC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(fcAUC))]);
    disp(['PC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pcAUC))]);
    disp(['WCS AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(wcsAUC))]);
    disp(['mvGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(gcAUC))]);
    disp(['pwGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pgcAUC))]);
    disp(['DLCM-GC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlAUC))]);
    disp(['LINUE-TE AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(linueAUC))]);
    disp(['DirectLiNGAM AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlgAUC))]);

    % save result
    fname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-result.mat'];
    save(fname, 'fcAUC','pcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlgAUC','linueAUC', 'fcROC','pcROC','wcsROC','gcROC','pgcROC','dlROC','dlgROC','linueROC');

    % show average ROC curve of DCM
    figure; 
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.9,0.5,0]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(pgcROC, N, [0.0,0.5,0.0]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
%    plotErrorROCcurve(dcmROC, N, [0.2,0.2,0.8]);
%    plotErrorROCcurve(rnnROC, N, [0.8,0.8,0.2]);
    plotErrorROCcurve(linueROC, N, [0.2,0.6,0.8]);
%    plotErrorROCcurve(nnnueROC, N, [0.8,0.2,0.8]);
    plotErrorROCcurve(dlgROC, N, [0.6,0.6,0.6]);
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '--', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '--', [0.0,0.5,0.0],0.5);
    plotAverageROCcurve(dlROC, N, '-', [0.2,0.2,0.2],1.2);
%    plotAverageROCcurve(dcmROC, N, '-', [0.2,0.2,0.8],0.5);
%    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5);
    plotAverageROCcurve(linueROC, N, '--', [0.2,0.5,0.7],0.5);
%    plotAverageROCcurve(nnnueROC, N, '--', [0.7,0.2,0.7],0.5);
    plotAverageROCcurve(dlgROC, N, '--', [0.6,0.6,0.6],0.5);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(i)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
