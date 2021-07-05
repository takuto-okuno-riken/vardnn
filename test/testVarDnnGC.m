
function testVarDnnGC
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
    % do training or load VARDNN network
    netFile = ['results/vardnn-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARDNN network
        net = initMvarDnnNetwork(si, exSignal, [], exControl);
        % training VARDNN network
        net = trainMvarDnnNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

        % recoverty training
        %[net, time] = recoveryTrainMvarDnnNetwork(si, exSignal, [], exControl, net, options);
        save(netFile, 'net');
    end
    
    % simulate VARDNN network with 1st frame & exogenous input signal
    [S, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, net);

    [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-DI)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);

    % show VARDNN-GC, VARDNN-DI
    figure; dlGC = plotMvarDnnGCI(si, exSignal, [], exControl, net, 0);
    figure; dlDI = plotMvarDnnDI(net, [], exControl, 0);
end

