
function testVarLstmGC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    lags = 1;
    % control is all positive input
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);

    % set training options
    maxEpochs = 1000;
    sigLen = size(si,2);
    miniBatchSize = ceil(sigLen / 2);

    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'SequenceLength','longest', ...
        'Shuffle','never', ...
        'GradientThreshold',5,...
        'Verbose',false);
%        'L2Regularization',0.005, ...
%        'Plots','training-progress', ...

    %% test pattern 1 
    % do training or load VARLSTM network
    netFile = ['results/mvarlstm-gc-test' num2str(nodeNum) '-' num2str(exNum) '-' num2str(lags) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARLSTM network
        net = initMvarLstmNetwork(si, exSignal, [], exControl);
        % training VARLSTM network
        net = trainMvarLstmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end

    % show VARLSTM-GC, VARLSTM-DI
    figure; dlGC = plotMvarLstmGCI(si, exSignal, [], exControl, net, 0);
    figure; dlDI = plotMvarLstmDI(net, [], exControl, 0);

    %% test pattern 2
    lags = 3;
    % do training or load VARLSTM network
    netFile = ['results/mvarlstm-gc-test' num2str(nodeNum) '-' num2str(exNum) '-' num2str(lags) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARLSTM network
        net = initMvarLstmNetwork(si, exSignal, [], exControl, lags);
        % training VARLSTM network
        net = trainMvarLstmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end

    % show VARLSTM-GC, VARLSTM-DI
    figure; dlGC = plotMvarLstmGCI(si, exSignal, [], exControl, net, 0);
    figure; dlDI = plotMvarLstmDI(net, [], exControl, 0, true);

    %% test pattern 3
    lags = 1;
    exNum = 4;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    exControl = logical(ones(nodeNum,exNum));
    si(3,2:end) = exSignal(2,1:sigLen-1);
    % do training or load VARLSTM network
    netFile = ['results/mvarlstm-gc-test' num2str(nodeNum) '-' num2str(exNum) '-' num2str(lags) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARLSTM network
        net = initMvarLstmNetwork(si, exSignal, [], exControl, lags);
        % training VARLSTM network
        net = trainMvarLstmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end

    % show VARLSTM-GC, VARLSTM-DI
    figure; dlGC = plotMvarLstmGCI(si, exSignal, [], exControl, net, 0, 0.05, true);
    figure; dlDI = plotMvarLstmDI(net, [], exControl, 0, true);
end

