
function testPvarDnnGC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
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
        'L2Regularization',0.05, ...
        'GradientThreshold',5,...
        'Verbose',false);
%            'Plots','training-progress');

    %% test pattern 1 
    % do training or load pairwise VAR DNN network
    lags = 1;
    pvardnnFile = ['results/pvardnn' num2str(lags) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pvardnnFile, 'file')
        load(pvardnnFile);
    else
        % init pairwise VAR DNN network
        netPvarDNN = initPvarDnnNetwork(si, exSignal, [], exControl, lags);
        % training pairwise VAR DNN network
        netPvarDNN = trainPvarDnnNetwork(si, exSignal, [], exControl, netPvarDNN, options);
        save(pvardnnFile, 'netPvarDNN');
    end
    
    % show pairwise DNN-GC
    figure; dnGC = plotPvarDnnGCI(si, exSignal, [], exControl, netPvarDNN, 0);
    % show pairwise DNN-DI
    figure; dnDI = plotPvarDnnDI(netPvarDNN, [], exControl);

    %% test pattern 2
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % do training or load pairwise VAR DNN network
    pvardnnFile = ['results/pvardnn' num2str(lags) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pvardnnFile, 'file')
        load(pvardnnFile);
    else
        % init pairwise VAR DNN network
        netPvarDNN = initPvarDnnNetwork(si, exSignal, [], exControl, lags);
        % training pairwise VAR DNN network
        netPvarDNN = trainPvarDnnNetwork(si, exSignal, [], exControl, netPvarDNN, options);
        save(pvardnnFile, 'netPvarDNN');
    end
    
    % show pairwise VAR DNN-GC
    figure; dnGC = plotPvarDnnGCI(si, exSignal, [], exControl, netPvarDNN, 0, 0.05, 1);
    % show pairwise VAR DNN-DI
    figure; dnDI = plotPvarDnnDI(netPvarDNN, [], exControl, 0, 1);

    %% test pattern 3
    exNum = 3;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    nodeControl = ones(nodeNum,nodeNum,lags);
    for i=1:nodeNum, nodeControl(i,i,2) = 0; end
    exControl = ones(nodeNum,exNum,lags);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    si(1,3:end) = exSignal(2,1:sigLen-2); % lag=2, this will be blocked by exControl
    exControl(1,2,2) = 0; % <= comment out and check control effect
    si(5,4:end) = si(1,1:sigLen-3); % lag=3, this will be blocked by nodeControl
    nodeControl(5,1,3) = 0; % <= comment out and check control effect

    % do training or load pairwise VAR DNN network
    pvardnnFile = ['results/pvardnn' num2str(lags) '-gc-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(pvardnnFile, 'file')
        load(pvardnnFile);
    else
        % init pairwise VAR DNN network
        netPvarDNN = initPvarDnnNetwork(si, exSignal, nodeControl, exControl, lags);
        % training pairwise VAR DNN network
        netPvarDNN = trainPvarDnnNetwork(si, exSignal, nodeControl, exControl, netPvarDNN, options);
        save(pvardnnFile, 'netPvarDNN');
    end
    
    % show pairwise VAR DNN-GC
    figure; dnGC = plotPvarDnnGCI(si, exSignal, nodeControl, exControl, netPvarDNN, 0, 0.05, 1);
    % show pairwise VAR DNN-DI
    figure; dnDI = plotPvarDnnDI(netPvarDNN, nodeControl, exControl, 0, 1);
end

