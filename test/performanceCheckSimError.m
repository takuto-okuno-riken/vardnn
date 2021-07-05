
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
%{
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
    bias = [0.5, 0, 0];
%}
%{
    load('data/marmoset-aneth-sample2-roi225.mat');
    siOrg = convert2SigmoidSignal(si);
    load('test/testTrain-rand500-uniform.mat');
    uuOrg = si;
    prefix = 'ms3';
    l2 = 0.005;
    imax = 1;
    inum = 450;
    weightFunc = @estimateInitWeightRoughHe;
    weightParam = [10];
    bias = [0.5, 0, 0];
%}
%%{
    load('data/marmoset-aneth-sample2-roi225.mat');
    siOrg = convert2SigmoidSignal(si);
    load('test/testTrain-rand500-uniform.mat');
    uuOrg = si;
    prefix = 'ms4';
    l2 = 0.005;
    imax = 1;
    inum = 450;
    weightFunc = [];
    weightParam = [];
    bias = 0;
%%}
    % do training and simulation and plot error graph
%%{
    exNum = 10;
    sigLen = 200;
    winLen = 200;
    for i=1:imax
        nodeNum = inum * i;

        si = siOrg(1:nodeNum,1:sigLen);
        exSignal = uuOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
        % control is all positive input
        exControl = logical(ones(nodeNum,exNum));

        % training and simulation
        checkingPattern(si, exSignal, exControl, winLen, i, prefix, l2, weightFunc, weightParam, bias);
    end
%%}
    % plot wisker-box graph
    Mae = []; R = []; FCcos = []; DLWcos = []; Tm = []; mTm = [];
    for i=1:imax
        nodeNum = inum * i;

        netFile = ['results/simerr/' prefix num2str(i) '_' num2str(nodeNum) 'x' num2str(winLen) '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
        load(netFile);
        Mae = [Mae; eachMae];
        R = [R; eachR];
        FCcos = [FCcos; eachFCcos];
        DLWcos = [DLWcos; eachDLWcos];
        mTm = [mTm, mean(allTime)];
        Tm = [Tm, allTime];
    end
    % box-and-whisker plot of MAE
    figure;
    boxplot(Mae);
    ylim([0 0.5]); title('box-and-whisker plot of MAE')
    % box-and-whisker plot of corr
    figure;
    boxplot(R);
    ylim([0 1]); title('box-and-whisker plot of corr')
    % box-and-whisker plot of cos similarity
    figure;
    boxplot(FCcos);
    ylim([0 1]); title('box-and-whisker plot of FC cos similarity');
    % box-and-whisker plot of cos similarity
    figure;
    boxplot(DLWcos);
    ylim([0 1]); title('box-and-whisker plot of VARDNN-DI cos similarity');
    % bar graph of simulation time of 100 signal length
    figure;
    bar(mTm);
end

function checkingPattern(si, exSignal, exControl, winLen, idx, prefix, l2, weightFunc, weightParam, bias)
    % traial number
    maxTrain = 8;
    maxWin = 1;
    % init
    nodeNum = size(si,1);
    exNum = size(exSignal,1);
    sigLen = size(si,2);
    allErr = []; allrS = []; allrSi = []; allTime = []; % all trained result
    allErr2 = []; 
    eachMae = []; eachR = []; eachFCcos = []; eachDLWcos = []; % each trained result
    ph = []; ch = [];

    for k = 1:maxTrain
        % do training or load VARDNN network
        netFile = ['results/simerr/' prefix num2str(idx) '-' num2str(k) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            % init VARDNN network
            netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl, 1, @reluLayer, weightFunc, weightParam, bias);

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

            % training VARDNN network
            netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
            % recover training 
            [netDLCM, time] = recoveryTrainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
            disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            %plotMvarDnnWeight(netDLCM);
            save(netFile, 'netDLCM');
        end
        
        allS2 = cell(maxWin,1);
        % simulate VARDNN network with 1st frame & exogenous input signal
        netFile = ['results/simerr/' prefix num2str(idx) '-' num2str(k) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(winLen) 'sim.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            allS = cell(maxWin,1);
            simTime = zeros(maxWin,1);
        end
        winSi=[]; winS=[]; winS2=[]; winFCcos=[]; winDLWcos=[];
        for i=1:maxWin
            st = floor(1 + (i-1)*(size(si,2)-winLen)/maxWin);
            en = st + winLen - 1;
            wSi = si(:,st:en);
            sExSignal = exSignal(:,st:en);
            % do simulation
            if isempty(allS{i})
                [S, time] = simulateMvarDnnNetwork(wSi, sExSignal, [], exControl, netDLCM);
                allS{i,1} = S; simTime(i) = time;
            else
                S = allS{i,1}; time = simTime(i);
                S2 = allS2{i,1};
            end
            if isempty(allS2{i})
                % simulate mPCVAR network with 1st frame & exogenous input signal
                netMVAR = initMpcvarNetwork(wSi, sExSignal, [], exControl, 2);
                [S2, time] = simulateMpcvarNetwork(wSi, sExSignal, [], exControl, netMVAR);
                allS2{i,1} = S2;
            else
                S2 = allS2{i,1};
            end
            % plot signals
            figure; plotTwoSignals(wSi(1:10,:), S(1:10,:)); title('DLCM simulation');
            figure; plotTwoSignals(wSi(1:10,:), S2(1:10,:)); title('PCVAR simulation');

            % keep for correlation
            winSi = [winSi; wSi]; winS = [winS; S]; winS2 = [winS2; S2];
            allrSi = [allrSi; wSi]; allrS = [allrS; S];
            % error of each window in one training
            [mae, maeerr, errs] = getTwoSignalsError(wSi, S);
            [mae2, maeerr, errs2] = getTwoSignalsError(wSi, S2);
            allErr  = [allErr; errs]; allErr2 = [allErr2; errs2];
            allTime = [allTime; time];
            % cosine similarity between FC and simulated FC
            FC = calcFunctionalConnectivity(wSi);
            smFC = calcFunctionalConnectivity(S); % simulated by DLCM
            smmpvFC = calcFunctionalConnectivity(S2); % simulated by PCVAR
            cs = getCosSimilarity(FC,smFC);
            cs2 = getCosSimilarity(FC,smmpvFC);
            winFCcos = [winFCcos; cs, cs2];
            % cosine similarity between original VARDNN-DI and simulated singal VARDNN-DI
            DLW = calcVarDnnConnectivity(wSi, sExSignal, [], exControl, '', k, i, winLen);
            smDLW = calcVarDnnConnectivity(S, sExSignal, [], exControl, 'sm', k, i, winLen);
            smmpvDLW = calcVarDnnConnectivity(S2, sExSignal, [], exControl, 'smmpv', k, i, winLen);
            cs = getCosSimilarity(DLW,smDLW);
            cs2 = getCosSimilarity(DLW,smmpvDLW);
            winDLWcos  = [winDLWcos; cs, cs2];
            disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
        end
        save(netFile, 'allS', 'simTime', 'allS2');
        % show error line graph
        Y = mean(abs(errs),1);
        Y2 = mean(abs(errs2),1);
        if isempty(ph)
            ph = figure;
        else
            figure(ph);
        end
        hold on;
        plot(Y, 'Color', [0.7,0.7,0.7]);
        plot(Y2, 'Color', [0.5,0.5,0.5]);
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
        figure; R = plotTwoSignalsCorrelation(winSi, winS) % show R result
        figure; R2 = plotTwoSignalsCorrelation(winSi, winS2) % show R result
        eachR = [eachR; R, R2];
        % error of each training
        [mae, maeerr, errs] = getTwoSignalsError(winSi, winS);
        [mae2, maeerr, errs] = getTwoSignalsError(winSi, winS2);
        eachMae = [eachMae; mae, mae2];
        % cos similarity of FC of each training
        eachFCcos = [eachFCcos; mean(winFCcos, 1)];
        % cos similarity of VARDNN-DI of each training
        eachDLWcos = [eachDLWcos; mean(winDLWcos, 1)];
    end
    % show mean all error line graph
    Y = mean(abs(allErr),1);
    Y2 = mean(abs(allErr2),1);
    figure(ph);
    hold on;
    plot(Y, 'Color', [0.2,0.2,1], 'LineWidth', 1);
    plot(Y2, 'Color', [1,0.2,0.2], 'LineWidth', 1);
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

    netFile = ['results/simerr/' prefix num2str(idx) '_' num2str(nodeNum) 'x' num2str(winLen)  '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
    save(netFile, 'allErr', 'allErr2', 'allrSi', 'allrS', 'allTime', 'eachMae', 'eachR', 'eachFCcos', 'eachDLWcos');
end

function mat = calcVarDnnConnectivity(X, exSignal, inControl, exControl, group, k, i, winLen)
    ROINUM = size(X,1);
    sigLen = size(X,2);
    netName = ['results/simerr/-dlcm-' num2str(ROINUM) 'x' num2str(winLen) '-' group num2str(k) '-' num2str(i) '.mat'];
    if exist(netName, 'file')
        f = load(netName);
        mat = f.mat;
    else
        netDLCM = initMvarDnnNetwork(X, exSignal, inControl, exControl);
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
        netDLCM = trainMvarDnnNetwork(X, exSignal, inControl, exControl, netDLCM, options);
        mat = calcMvarDnnDI(netDLCM, inControl, exControl);
        save(netName, 'netDLCM','mat');
    end
end
