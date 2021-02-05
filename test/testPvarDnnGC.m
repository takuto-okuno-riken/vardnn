
function testPvarDnnGC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 3;
    nodeNum = 6;
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
        'L2Regularization',0.01, ...
        'GradientThreshold',1,...
        'Verbose',false);
%            'Plots','training-progress');

    %% test pattern 1 
    % do training or load Pairwised VAR DNN network
    pvardnnFile = ['results/pvardnn-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pvardnnFile, 'file')
        load(pvardnnFile);
    else
        % init Pairwised VAR DNN network
        netPvarDNN = initPvarDnnNetwork(si, exSignal, [], exControl, lags);
        % training Pairwised VAR DNN network
        netPvarDNN = trainPvarDnnNetwork(si, exSignal, [], exControl, netPvarDNN, options);
        save(pvardnnFile, 'netPvarDNN');
    end
    
    % show Pairwised DNN-GC
    figure; dnGC = plotPvarDnnGCI(si, exSignal, [], exControl, netPvarDNN, 0);
    % show Pairwised DNN-EC
    figure; dnEC = plotPvarDnnEC(netPvarDNN, [], exControl);

    %% test pattern 2
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % do training or load Pairwised VAR DNN network
    pvardnnFile = ['results/pvardnn-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pvardnnFile, 'file')
        load(pvardnnFile);
    else
        % init Pairwised VAR DNN network
        netPvarDNN = initPvarDnnNetwork(si, exSignal, [], exControl, lags);
        % training Pairwised VAR DNN network
        netPvarDNN = trainPvarDnnNetwork(si, exSignal, [], exControl, netPvarDNN, options);
        save(pvardnnFile, 'netPvarDNN');
    end
    
    % show Pairwised VAR DNN-GC
    figure; dnGC = plotPvarDnnGCI(si, exSignal, [], exControl, netPvarDNN, 0, 0.05, 1);
    % show Pairwised VAR DNN-EC
    figure; dnEC = plotPvarDnnEC(netPvarDNN, [], exControl, 0, 1);
end

