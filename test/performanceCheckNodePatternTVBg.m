function performanceCheckNodePatternTVBg
    num_scan = 55;
    if num_scan == 55  % oh's mouse 98 node. density 0.15. weight add.
        Gvals = {'0.25', '0.5', '1.0', '1.5', '2.0'};
        node_nums = [16,98];
        GTths = [1, 1];
        pat = [1, 6];
    end
    % test sparse and full density
    hz = 64;
    N = 8;
    for k=1:length(Gvals)
        for i=1:length(node_nums)
            checkingPattern(Gvals{k}, node_nums(i), num_scan, hz, GTths(i), N, pat(i));
        end
    end
end

function checkingPattern(Gval, node_num, num_scan, hz, GTth, N, i)
    % init
    fcAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    fcROC = cell(N,2);
    gcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);

    origf = figure;
    origSigf = figure;

    lag = 3;

    for k=1:N
        tvbFile = ['data/rww55g/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan' Gval 'g-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
        load(tvbFile);
        density = length(find(weights>GTth)) / (node_num * (node_num-1));
        nodeNum = size(si,1);
        sigLen = size(si,2);

        % show original connection
        figure(origf); plotEC(weights, 'Ground Truth', 1);
        figure(origSigf); plot(t, si);

        % show result of FC
        FC = calcFunctionalConnectivity(si);
        [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = calcROCcurve(FC, weights, 100, 1, GTth);

        % show result of granger causality index (mvGC)
        gcI = calcMultivariateGCI_(si,[],[],[],lag);
        [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = calcROCcurve(gcI, weights, 100, 1, GTth);
%%{
        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(si);
        [uu, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            
        % calcurate and show DLCM-GC
        dlGC = [];
        exControl = eye(nodeNum, nodeNum);
        netFile = ['results/rww55g/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-' Gval 'g-idx' num2str(i) '-' num2str(k) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
            if exist('inSignal','var'), exSignal=inSignal; end % for compatibility
        else
            % train VARDNN    
            Y = si;
            exSignal = uu;
            % layer parameters
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

            disp('start training');
            netDLCM = trainMvarDnnNetwork(Y, exSignal, [], exControl, netDLCM, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);

            save(netFile, 'netDLCM', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
        end
        if isempty(dlGC)
            % show DLCM-GC
            dlGC = calcMvarDnnGCI(Y, exSignal, [], exControl, netDLCM);
            save(netFile, 'netDLCM', 'Y', 'exSignal', 'si', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2', 'dlGC');
        end
        
        % calc ROC curve
        [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = calcROCcurve(dlGC, weights, 100, 1, GTth);

        % show result of VARDNN weight causality index (DLCM-wci) as DLCM-EC
        dlwGC = calcMvarDnnEC(netDLCM, [], exControl);
       [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = calcROCcurve(dlwGC, weights, 100, 1, GTth);
%%}
    end
    % show result AUC
    disp(['FC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(fcAUC))]);
    disp(['mvGC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(gcAUC))]);
    disp(['DLCM-GC AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlAUC))]);
    disp(['DLCM-WCI AUC (' num2str(i) ', node=' num2str(node_num) ', density=' num2str(density) ') : ' num2str(mean(dlwAUC))]);

    % save result
    fname = ['results/rww55g/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan' Gval 'g-pat' num2str(i) '-' num2str(hz) 'hz-result.mat'];
    save(fname, 'fcAUC','gcAUC','dlAUC','dlwAUC','fcROC','gcROC','dlROC','dlwROC');

    % show average ROC curve of DCM
    figure; 
    hold on;
%    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % TODO:
%    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.2);
    plotAverageROCcurve(dlROC, N, '--', [0.2,0.2,0.2],0.7); % TODO:
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve G=' Gval ', node=' num2str(node_num)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
