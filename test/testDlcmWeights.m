
function testDlcmWeights
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
    netFile = ['results/dlcm-w-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
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

    % show VARDNN representative weights ------------------------------------
    sampling = eye(nodeNum);
    sampling2 = zeros(nodeNum);
    repW = zeros(nodeNum,nodeNum);
    repW2 = zeros(nodeNum,nodeNum);
    for i=1:nodeNum
        % predict representative weights
        repW(i,:) = predict(netDLCM.nodeNetwork{i}, sampling);
        repW2(i,:) = predict(netDLCM.nodeNetwork{i}, sampling2);
    end
    %figure; imagesc(repW); daspect([1 1 1]);
    %figure; imagesc(repW2); daspect([1 1 1]);

    r2 = repW - repW2;
    figure; imagesc(r2); daspect([1 1 1]); colorbar;

    % show VARDNN representative weights 2 ----------------------------------
    allone = ones(nodeNum, 1);
    sampling3 = ones(nodeNum, nodeNum) - eye(nodeNum);
    repW3 = zeros(nodeNum,nodeNum);
    for i=1:nodeNum
        % predict representative weights
        all = predict(netDLCM.nodeNetwork{i}, allone);
        R = predict(netDLCM.nodeNetwork{i}, sampling3);
        repW3(i,:) = all - R;
    end
    figure; imagesc(repW3); daspect([1 1 1]); colorbar;
    
    %
    r4 = (r2 + repW3) / 2;
    figure; imagesc(r4); daspect([1 1 1]); colorbar;
    
    % VARDNN weight causal index as VARDNN-DI
    wci = plotMvarDnnEC(netDLCM, [], [], 0);
end

