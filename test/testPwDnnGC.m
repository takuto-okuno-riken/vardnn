
function testPwDnnGC
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
    % do training or load Pairwised DNN-GC's network
    pwdnnFile = ['results/pw-dnn-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pwdnnFile, 'file')
        load(pwdnnFile);
    else
        % init Pairwised DNN-GC's network
        netPwDNN = initPwDnnGCNetwork(si, exSignal, [], exControl, lags);
        % training Pairwised DNN-GC's network
        netPwDNN = trainPwDnnGCNetwork(si, exSignal, [], exControl, netPwDNN, options);
        save(pwdnnFile, 'netPwDNN');
    end
    
    % show Pairwised DNN-GC
    figure; dnGC = plotPwDnnGCI(si, exSignal, [], exControl, netPwDNN, 0);


    %% test pattern 2
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % do training or load Pairwised DNN-GC's network
    pwdnnFile = ['results/pw-dnn-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pwdnnFile, 'file')
        load(pwdnnFile);
    else
        % init Pairwised DNN-GC's network
        netPwDNN = initPwDnnGCNetwork(si, exSignal, [], exControl, lags);
        % training Pairwised DNN-GC's network
        netPwDNN = trainPwDnnGCNetwork(si, exSignal, [], exControl, netPwDNN, options);
        save(pwdnnFile, 'netPwDNN');
    end
    
    % show Pairwised DNN-GC
    figure; dnGC = plotPwDnnGCI(si, exSignal, [], exControl, netPwDNN, 0, 0.05, 1);
end

