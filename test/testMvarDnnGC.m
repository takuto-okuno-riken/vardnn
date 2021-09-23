
function testMvarDnnGC
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
    % do training or load multivariate VAR DNN network
    for i=1:5
        mvardnnFile = ['results/mvardnn' num2str(i) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(mvardnnFile, 'file')
            load(mvardnnFile);
        else
            % init multivariate VAR DNN network
            net = initMvarDnnNetwork(si, exSignal, [], exControl, i);
            % training multivariate VAR DNN network
            net = trainMvarDnnNetwork(si, exSignal, [], exControl, net, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(net);
            disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(mvardnnFile, 'net');
        end

        % show mVAR-DNN-GC
        figure; mvdnnGC = plotMvarDnnGCI(si, exSignal, [], exControl, net, 0);
        figure; mvdnnDI = plotMvarDnnDI(net, [], exControl, 0);
        figure; mvdnnMIV = plotMvarDnnMIV(si, exSignal, [], exControl, net, 0);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], exControl, i, 0);
    end

    %% test pattern 2 -- exogenous signals
    % do training or load multivariate VAR DNN network
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    for i=1:5
        mvardnnFile = ['results/mvardnn' num2str(i) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(mvardnnFile, 'file')
            load(mvardnnFile);
        else
            % init multivariate VAR DNN network
            net = initMvarDnnNetwork(si, exSignal, [], exControl, i);
            % training multivariate VAR DNN network
            net = trainMvarDnnNetwork(si, exSignal, [], exControl, net, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(net);
            disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(mvardnnFile, 'net');
        end

        % show mVAR-DNN-GC
        figure; mvdnnGC = plotMvarDnnGCI(si, exSignal, [], exControl, net, 0, 0, 1);
        figure; mvdnnDI = plotMvarDnnDI(net, [], exControl, 0, 1);
        figure; mvdnnMIV = plotMvarDnnMIV(si, exSignal, [], exControl, net, 0, 1);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], exControl, i, 0, 0, 1);
    end
    
    %% test pattern 3 -- no activation function
    % do training or load multivariate VAR DNN network
    exNum = 3;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    for i=1:5
        mvardnnFile = ['results/mvardnn' num2str(i) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(mvardnnFile, 'file')
            load(mvardnnFile);
        else
            % init multivariate VAR DNN network. no activation function
            net = initMvarDnnNetwork(si, exSignal, [], exControl, i, []);
            % training multivariate VAR DNN network
            net = trainMvarDnnNetwork(si, exSignal, [], exControl, net, options);
            [time, loss, rsme] = getMvarDnnTrainingResult(net);
            disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(mvardnnFile, 'net');
        end

        % show mVAR-DNN-GC
        figure; mvdnnGC = plotMvarDnnGCI(si, exSignal, [], exControl, net, 0, 0, 1);
        figure; mvdnnDI = plotMvarDnnDI(net, [], exControl, 0, 1);
        figure; mvdnnMIV = plotMvarDnnMIV(si, exSignal, [], exControl, net, 0, 1);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], exControl, i, 0, 0, 1);
    end
end

