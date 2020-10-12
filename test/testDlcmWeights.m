
function testDlcmWeights
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    inputNum = 0;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    inSignal = [];
    % control is all positive input
    inControl = [];
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
    % do training or load DLCM network
    dlcmFile = ['results/dlcm-w-test' num2str(nodeNum) '-' num2str(inputNum) '.mat'];
    if exist(dlcmFile, 'file')
        load(dlcmFile);
    else
        % init DLCM network
        netDLCM = initDlcmNetwork(si, inSignal, [], inControl);
        % training DLCM network
        netDLCM = trainDlcmNetwork(si, inSignal, [], inControl, netDLCM, options);
        [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

        % recoverty training
        %[netDLCM, time] = recoveryTrainDlcmNetwork(si, inSignal, [], inControl, netDLCM, options);
        save(dlcmFile, 'netDLCM');
    end

    % show DLCM representative weights ------------------------------------
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

    % show DLCM representative weights 2 ----------------------------------
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
    
    % DLCM weight causal index as DLCM-EC
    wci = plotDlcmEC(netDLCM, [], [], 0);
end

