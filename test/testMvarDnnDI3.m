
function testMvarDnnDI3
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    lags = 1;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    nodeControl = [];
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1); % lag=1
    si(4,3:end) = si(6,1:sigLen-2); % lag=2

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
    % do training or load multivariate VAR DNN network
    netFile = ['results/mvardnn' num2str(lags) '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init multivariate VAR DNN network
        net = initMvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR DNN network
        net = trainMvarDnnNetwork(si, exSignal, nodeControl, exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end

    % show mVAR-DNN-DI3
    DI3 = calcMvarDnnDI3(net, nodeControl, exControl, 0);
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, nodeControl, exControl, net, 0);
    figure; plotDirectedFC3(DI3,'mVARDNN-DI',0);
    figure; plotDirectedFC3(MIV3,'mVARDNN-MIV',0);
    figure; plotDirectedFC3(MAIV3,'mVARDNN-MAIV',0);

    %% test pattern 2 -- exogenous signals
    % do training or load multivariate VAR DNN network
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum,lags);
    si(3,2:end) = exSignal(1,1:sigLen-1); % lag=1
    si(5,4:end) = exSignal(2,1:sigLen-3); % lag=3
    netFile = ['results/mvardnn' num2str(lags) '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init multivariate VAR DNN network
        net = initMvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR DNN network
        net = trainMvarDnnNetwork(si, exSignal, nodeControl, exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end

    % show mVAR-DNN-GC
    DI3 = calcMvarDnnDI3(net, nodeControl, exControl, 1);
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, nodeControl, exControl, net, 1);
    figure; plotDirectedFC3(DI3,'mVARDNN-DI',0);
    figure; plotDirectedFC3(MIV3,'mVARDNN-MIV',0);
    figure; plotDirectedFC3(MAIV3,'mVARDNN-MAIV',0);

    %% test pattern 3 -- input value selection test
    % do training or load multivariate VAR DNN network
    lags = 10;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum,lags);
    si(3,4:end) = exSignal(1,1:sigLen-3); % lag=3
    si(5,8:end) = exSignal(2,1:sigLen-7); % lag=7
    netFile = ['results/mvardnn' num2str(lags) '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init multivariate VAR DNN network
        net = initMvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR DNN network
        net = trainMvarDnnNetwork(si, exSignal, nodeControl, exControl, net, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net');
    end
    % input value selection by VARDNN-DI
    ivsName = 'ivs1';
    rate = 8 / (10*10*lags);
    DI3 = calcMvarDnnDI3(net, [], [], 1);
    control3 = inputValueSelectionFromDFC3(DI3, rate);
    nodeControl = control3(:,1:nodeNum,:);
    exControl = control3(:,nodeNum+1:nodeNum+exNum,:);
    netFile = ['results/mvardnn' num2str(lags) '-' ivsName '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init multivariate VAR DNN network
        net2 = initMvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR DNN network
        net2 = trainMvarDnnNetwork(si, exSignal, nodeControl, exControl, net2, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net2);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net2');
    end
    DI3 = calcMvarDnnDI3(net2, nodeControl, exControl, 1);    
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, nodeControl, exControl, net2, 1);
    figure; plotDirectedFC3(DI3,[ivsName ' mVARDNN-DI'],0);
    figure; plotDirectedFC3(MIV3,[ivsName ' mVARDNN-MIV'],0);
    figure; plotDirectedFC3(MAIV3,[ivsName ' mVARDNN-MAIV'],0);

    % input value selection by VARDNN-MIV
    ivsName = 'ivs2';
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, [], [], net, 1);
    control3 = inputValueSelectionFromDFC3(MIV3, rate);
    nodeControl = control3(:,1:nodeNum,:);
    exControl = control3(:,nodeNum+1:nodeNum+exNum,:);
    netFile = ['results/mvardnn' num2str(lags) '-' ivsName '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init multivariate VAR DNN network
        net2 = initMvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR DNN network
        net2 = trainMvarDnnNetwork(si, exSignal, nodeControl, exControl, net2, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net2);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net2');
    end
    DI3 = calcMvarDnnDI3(net2, nodeControl, exControl, 1);    
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, nodeControl, exControl, net2, 1);
    figure; plotDirectedFC3(DI3,[ivsName ' mVARDNN-DI'],0);
    figure; plotDirectedFC3(MIV3,[ivsName ' mVARDNN-MIV'],0);
    figure; plotDirectedFC3(MAIV3,[ivsName ' mVARDNN-MAIV'],0);

    % input value selection by VARDNN-DI
    ivsName = 'ivs2';
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, [], [], net, 1);
    control3 = inputValueSelectionFromDFC3(MAIV3, rate);
    nodeControl = control3(:,1:nodeNum,:);
    exControl = control3(:,nodeNum+1:nodeNum+exNum,:);
    netFile = ['results/mvardnn' num2str(lags) '-' ivsName '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init multivariate VAR DNN network
        net2 = initMvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training multivariate VAR DNN network
        net2 = trainMvarDnnNetwork(si, exSignal, nodeControl, exControl, net2, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(net2);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        save(netFile, 'net2');
    end
    DI3 = calcMvarDnnDI3(net2, nodeControl, exControl, 1);    
    [MIV3,MAIV3] = calcMvarDnnMIV3(si, exSignal, nodeControl, exControl, net2, 1);
    figure; plotDirectedFC3(DI3,[ivsName ' mVARDNN-DI'],0);
    figure; plotDirectedFC3(MIV3,[ivsName ' mVARDNN-MIV'],0);
    figure; plotDirectedFC3(MAIV3,[ivsName ' mVARDNN-MAIV'],0);
end



