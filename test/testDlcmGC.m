
function testDlcmGC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);

    % set training options
    maxEpochs = 1000;
    sigLen = size(si,2);
    miniBatchSize = ceil(sigLen / 3);

    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'Shuffle','every-epoch', ...
        'L2Regularization',0.05, ...
        'GradientThreshold',5,...
        'Verbose',false);
%            'Plots','training-progress');

    %% test pattern 1 
    % do training or load DLCM network
    netFile = ['results/dlcm-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARDNN network
        netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
        % training VARDNN network
        netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

        % recoverty training
        %[netDLCM, time] = recoveryTrainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
        save(netFile, 'netDLCM');
    end
    
    % simulate DLCM network with 1st frame & exogenous input signal
    [S, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, netDLCM);

    [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);

    % show DLCM-GC, DLCM-EC
    figure; dlGC = plotMvarDnnGCI(si, exSignal, [], exControl, netDLCM, 0);
    figure; dlEC = plotMvarDnnEC(netDLCM, [], exControl, 0);
end

