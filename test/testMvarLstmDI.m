
function testMvarLstmDI
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 1;
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
        'SequenceLength','longest', ...
        'Shuffle','never', ...
        'GradientThreshold',5,...
        'Verbose',false);
%            'Plots','training-progress');

    %% test pattern 1 
    netFile = ['results/mvarlstm' num2str(lags) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        net = initMvarLstmNetwork(si, exSignal, [], exControl, lags);
        % training multivariate VAR LSTM network
        net = trainMvarLstmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end
    % show multivaliate VARLSTM-DI
    figure; mDI = plotMvarLstmDI(net, [], exControl, 0);
    % show multivaliate VARLSTM-GC
    figure; mGC = plotMvarLstmGCI(si, exSignal, [], exControl, net, 0);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0);

    
    %% test pattern 2
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % init VARLSTM network
    netFile = ['results/mvarlstm' num2str(lags) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        net = initMvarLstmNetwork(si, exSignal, [], exControl, lags);
        % training multivariate VAR LSTM network
        net = trainMvarLstmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end
    % show multivaliate VARLSTM-DI
    figure; mDI = plotMvarLstmDI(net, [], exControl, 0, 1);
    % show multivaliate VARLSTM-GC
    figure; mGC = plotMvarLstmGCI(si, exSignal, [], exControl, net, 0, 0.05, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0, 0.05, 1);

    %% test pattern 3
    exNum = 4;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    nodeControl = ones(nodeNum,nodeNum,lags);
    for k=1:nodeNum, nodeControl(k,k,2)=0; end
    exControl = ones(nodeNum,exNum,lags);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    si(5,3:end) = si(1,1:sigLen-2); % lag=2, this will be blocked by nodeControl
    nodeControl(5,1,2) = 0; % <= comment out and check control effect
    si(7,3:end) = exSignal(2,1:sigLen-2); % lag=2, this will be blocked by exControl
    exControl(7,2,2) = 0; % <= comment out and check control effect

    % init VARLSTM network
    netFile = ['results/mvarlstm' num2str(lags) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        net = initMvarLstmNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR LSTM network
        net = trainMvarLstmNetwork(si, exSignal, nodeControl, exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end
    % show multivaliate VARLSTM-DI
    figure; mDI = plotMvarLstmDI(net, nodeControl, exControl, 0, 1);
    % show multivaliate VARLSTM-GC
    figure; mGC = plotMvarLstmGCI(si, exSignal, nodeControl, exControl, net, 0, 0.05, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, nodeControl, exControl, lags, 0, 0.05, 1);
end

