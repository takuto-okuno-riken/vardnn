
function performanceCheckSimError
    % load signals
%{
    load('test/testTrain-rand500-uniform.mat');
    prefix = '';
    l2 = 0.0005;
    siOrg = si;
    uuOrg = si;
    imax = 6;
    inum = 30;
    weightFunc = [];
    weightParam = [];
    bias = 0;
%}
%%{
    load('data/marmoset-aneth-sample2-roi225.mat');
    siOrg = bold2dnnSignal(si, 0.2);
    load('test/testTrain-rand500-uniform.mat');
    uuOrg = si;
    prefix = 'ms2';
    l2 = 0.005;
    imax = 1;
    inum = 450;
    weightFunc = @estimateInitWeightRoughHe;
    weightParam = [10];
    bias = 0.5;
%%}
    % do training and simulation and plot error graph
%%{
    for i=1:imax
        nodeNum = inum * i;
        exNum = 10;
        sigLen = 200;
        winLen = 100;

        si = siOrg(1:nodeNum,1:sigLen);
        exSignal = uuOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
        % control is all positive input
        exControl = logical(ones(nodeNum,exNum));

        % training and simulation
        checkingPattern(si, exSignal, exControl, winLen, i, prefix, l2, weightFunc, weightParam, bias);
    end
%%}
    % plot wisker-box graph
    Mae = []; R = []; FCcos = []; GCcos = []; Tm = []; mTm = [];
    for i=1:imax
        nodeNum = inum * i;
        exNum = 10;
        sigLen = 200;
        winLen = 100;

        dlcmFile = ['results/sim-err' prefix num2str(i) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
        load(dlcmFile);
        Mae = [Mae, eachMae];
        R = [R, eachR];
        FCcos = [FCcos, eachFCcos];
        GCcos = [GCcos, eachGCcos];
        mTm = [mTm, mean(allTime)];
        Tm = [Tm, allTime];
    end
    % box-and-whisker plot of MAE
    figure;
    boxplot(Mae);
    ylim([0 0.5]);
    % box-and-whisker plot of corr
    figure;
    boxplot(R);
    ylim([0 1]);
    % box-and-whisker plot of cos similarity
    figure;
    boxplot(FCcos);
    ylim([0 1]);
    % box-and-whisker plot of cos similarity
    figure;
    boxplot(GCcos);
    ylim([0 1]);
    % bar graph of simulation time of 100 signal length
    figure;
    bar(mTm);
end

function checkingPattern(si, exSignal, exControl, winLen, idx, prefix, l2, weightFunc, weightParam, bias)
    % traial number
    maxTrain = 8;
    maxWin = 10;
    % init
    nodeNum = size(si,1);
    exNum = size(exSignal,1);
    sigLen = size(si,2);
    allErr = []; allrS = []; allrSi = []; allTime = []; % all trained result
    eachMae = []; eachR = []; eachFCcos = []; eachGCcos = []; % each trained result
    ph = []; ch = [];

    for k = 1:maxTrain
        % do training or load DLCM network
        dlcmFile = ['results/net-sim-err' prefix num2str(idx) '-' num2str(k) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            % init DLCM network
            netDLCM = initDlcmNetwork(si, exSignal, [], exControl, weightFunc, weightParam, bias);

            % set training options
            maxEpochs = 1000;
            sigLen = size(si,2);
            miniBatchSize = ceil(sigLen / 3);

            options = trainingOptions('adam', ...
                'ExecutionEnvironment','cpu', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Shuffle','every-epoch', ...
                'GradientThreshold',5,...
                'L2Regularization',l2, ...
                'Verbose',false);
        %            'Plots','training-progress');

            % training DLCM network
            netDLCM = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
            % recover training 
            [netDLCM, time] = recoveryTrainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            %plotDlcmWeight(netDLCM);
            save(dlcmFile, 'netDLCM');
        end

        % simulate DLCM network with 1st frame & exogenous input signal
        dlcmFile = ['results/net-sim-err' prefix num2str(idx) '-' num2str(k) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(sigLen) 'sim.mat'];
        if exist(dlcmFile, 'file')
            load(dlcmFile);
        else
            allS = cell(maxWin,1);
            simTime = zeros(maxWin,1);
        end
        winSi=[]; winS=[]; winFCcos=[]; winGCcos=[];
        for i=1:maxWin
            st = floor(1 + (i-1)*(size(si,2)-winLen)/maxWin);
            en = st + winLen - 1;
            wSi = si(:,st:en);
            sExSignal = exSignal(:,st:en);
            % do simulation
            if isempty(allS{i})
                [S, time] = simulateDlcmNetwork(wSi, sExSignal, [], exControl, netDLCM);
                allS{i,1} = S; simTime(i) = time;
            else
                S = allS{i,1}; time = simTime(i);
            end

            % keep for correlation
            winSi = [winSi; wSi]; winS = [winS; S];
            allrSi = [allrSi; wSi]; allrS = [allrS; S];
            % error of each window in one training
            [mae, maeerr, errs] = getTwoSignalsError(wSi, S);
            allErr = [allErr; errs];
            allTime = [allTime; time];
            % cosine similarity between FC and simulated FC
            FC = calcFunctionalConnectivity(wSi);
            sFC = calcFunctionalConnectivity(S);
            cs = getCosSimilarity(FC,sFC);
            winFCcos = [winFCcos; cs];
            % cosine similarity between GC and simulated GC
            gcI  = calcPairwiseGCI(wSi, sExSignal, [], exControl, 3);
            sgcI = calcPairwiseGCI(S, sExSignal, [], exControl, 3);
            nidx = find(isnan(gcI)); gcI(nidx) = 0; % remove NaN
            nidx = find(isnan(sgcI)); sgcI(nidx) = 0; % remove NaN
            cs = getCosSimilarity(gcI,sgcI);
            winGCcos = [winGCcos; cs];
            disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
        end
        save(dlcmFile, 'allS', 'simTime');
        % show error line graph
        Y = mean(abs(errs),1);
        if isempty(ph)
            ph = figure;
        else
            figure(ph);
        end
        hold on;
        plot(Y, 'Color', [0.7,0.7,0.7]);
        ylim([0 0.5]);
        hold off;

        % show correlation line graph
        Y = zeros(1,winLen);
        for i=1:winLen
            Y(i) = corr2(winSi(:,i),winS(:,i));
        end
        if isempty(ch)
            ch = figure;
        else
            figure(ch);
        end
        hold on;
        plot(Y, 'Color', [0.7,0.7,0.7]);
        ylim([0 1]);
        hold off;

        % show correlation graph of each training
        figure;
        R = plotTwoSignalsCorrelation(winSi, winS) % show R result
        eachR = [eachR; R];
        % error of each training
        [mae, maeerr, errs] = getTwoSignalsError(winSi, winS);
        eachMae = [eachMae; mae];
        % cos similarity of FC of each training
        eachFCcos = [eachFCcos; mean(winFCcos)];
        % cos similarity of GC of each training
        eachGCcos = [eachGCcos; mean(winGCcos)];
    end
    % show mean all error line graph
    Y = mean(abs(allErr),1);
    figure(ph);
    hold on;
    plot(Y, 'Color', [0.2,0.2,1], 'LineWidth', 1);
    hold off;
    
    % show all correlation line graph
    Y = zeros(1,winLen);
    for i=1:winLen
        Y(i) = corr2(allrSi(:,i),allrS(:,i));
    end
    figure(ch);
    hold on;
    plot(Y, 'Color', [0.2,0.2,1], 'LineWidth', 1);
    hold off;
    drawnow;

    dlcmFile = ['results/sim-err' prefix num2str(idx) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
    save(dlcmFile, 'allErr', 'allrSi', 'allrS', 'allTime', 'eachMae', 'eachR', 'eachFCcos', 'eachGCcos');
end

