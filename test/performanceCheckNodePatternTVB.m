function performanceCheckNodePatternTVB
%    node_nums = [16, 16, 32, 32, 76];
    %node_nums = [16, 16, 16, 32, 32, 32];
    %node_nums = [8, 8, 8, 8, 8, 8, 8, 8];
%    node_nums = [16, 16, 16, 16, 16, 16];
%    node_nums = [32, 32, 32, 32, 32, 32];
%    node_nums = [76];
%    node_nums = [16, 16, 32, 32, 76];
    node_nums = [16, 32, 48, 64, 80, 96];
%    Gths = [0.2, 0.2, 0.2, 0.2, 0.2];
%    Gths = [0, 0, 0, 0, 0, 0];
    Gths = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2];
    %Gths = [0, 0, 0, 0, 0.2, 0.2, 0.2, 0.2];
%    Gths = [0.2];
    num_scan = 45;
    hz = 64;
    N = 8;

    for i=1:length(node_nums)
        checkingPattern(node_nums(i), num_scan, hz, Gths(i), N, i);
    end
end

function checkingPattern(node_num, num_scan, hz, Gth, N, i)
    % init
    fcAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    rnnAUC = zeros(1,N);
    linueAUC = zeros(1,N);
    nnnueAUC = zeros(1,N);
    fcROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    rnnROC = cell(N,2);
    linueROC = cell(N,2);
    nnnueROC = cell(N,2);
    fcRf = figure;
    gcRf = figure;
    pgcRf = figure;
    dlRf = figure;
%    rnnRf = figure;
    linueRf = figure;
%    nnnueRf = figure;
    origf = figure;
    origSigf = figure;

    rnnTrial = 8;
    lag = 3;

    for k=1:N
        tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
        load(tvbFile);

        [si, sig, m, maxsi, minsi] = convert2SigmoidSignal(si);
        [uu, sig2, m2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            
        % show original connection
        figure(origf); plotDcmEC(weights);
        figure(origSigf); plot(t, si);

        % show original signal FC
        FC = calcFunctionalConnectivity(si);
        figure(fcRf); hold on; [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = plotROCcurve(FC, weights, 100, 1, Gth); hold off;
        title(['ROC curve of FC (pat=' num2str(i) ')']);

        % show original signal granger causality index (mvGC)
        gcI = calcMultivariateGCI(si,lag);
        figure(gcRf); hold on; [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
        title(['ROC curve of mvGC (pat=' num2str(i) ')']);

        % show original signal granger causality index (pwGC)
        gcI = calcPairwiseGCI(si,lag);
        figure(pgcRf); hold on; [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
        title(['ROC curve of pwGC (pat=' num2str(i) ')']);
%%{
        % calcurate and show DLCM-GC
        nodeNum = size(si,1);
        sigLen = size(si,2);
        inControl = eye(nodeNum, nodeNum);
        dlcmFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            % train DLCM    
            %[Y, sig, m, maxsi, minsi] = convert2SigmoidSignal(si);
            %[inSignal, sig2, m2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            Y = si;
            inSignal = uu;
            % layer parameters
            weightFunc = @estimateInitWeightRoughHe;
            weightParam = [10];
            bias = 0.5;
            netDLCM = initDlcmNetwork(Y, inSignal, inControl); % weightFunc, weightParam, bias);
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
            netDLCM = trainDlcmNetwork(Y, inSignal, inControl, netDLCM, options);
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);

            % recoverty training
            %[netDLCM, time] = recoveryTrainDlcmNetwork(Y, inSignal, inControl, netDLCM, options);
            save(dlcmFile, 'netDLCM', 'Y', 'inSignal', 'Y', 'sig', 'm', 'maxsi', 'minsi', 'sig2', 'm2', 'maxsi2', 'minsi2');
        end
        % show DLCM-GC
        dlGC = calcDlcmGCI(Y, inSignal, inControl, netDLCM);
        
        % calc ROC curve
        figure(dlRf); hold on; [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = plotROCcurve(dlGC, weights); hold off;
        title(['ROC curve of DLCM-GC (pat=' num2str(i) ')']);
%%}
%{
        % linue TE result
        linueFile = ['results/tvb-pat-linue/linue_MultivAnalysis_tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '-' num2str(lag) '.mat'];
        load(linueFile);
        A = outputToStore.reshapedMtx.';

        % show ROC curve of TE(LIN UE)
        figure(linueRf); hold on; [linueROC{k,1}, linueROC{k,2}, linueAUC(k)] = plotROCcurve(A, weights, 100, 1, Gth); hold off;        
        title(['ROC curve of LINUE-TE (pat=' num2str(i) ')']);
%}
    end
    % show result AUC
    disp(['FC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(fcAUC))]);
    disp(['mvGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(gcAUC))]);
    disp(['pwGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(pgcAUC))]);
    disp(['DLCM-GC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlAUC))]);
    disp(['LINUE-TE AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(linueAUC))]);

    % save result
    fname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-result.mat'];
    save(fname, 'fcAUC','gcAUC','pgcAUC','dlAUC','rnnAUC','linueAUC','nnnueAUC', 'fcROC','gcROC','pgcROC','dlROC','rnnROC','linueROC','nnnueROC');
end
