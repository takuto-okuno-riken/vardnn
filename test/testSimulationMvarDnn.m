
function testSimulationMVAR
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 4;
    sigLen = 100;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = logical(ones(nodeNum,exNum));
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
        'GradientThreshold',5,...
        'L2Regularization',0.05, ...
        'Verbose',false);

    %% test pattern 1 
    % do training or load mVar DNN network
    lags = 1;
    mvardnnFile = ['results/mvardnn-sim-test8-4-' num2str(lags) '.mat'];
    if exist(mvardnnFile, 'file')
        load(mvardnnFile);
    else
        % init mVAR DNN network
        net = initDlcmNetwork(si, exSignal, [], exControl, lags);
        % training mVAR DNN network
        net = trainDlcmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getDlcmTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

        save(mvardnnFile, 'net');
    end
    
    % simulate mVAR DNN (DLCM) network with 1st frame & exogenous input signal
    [S, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, net);
    figure; [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    [S2, time] = simulateDlcmNetwork(si, exSignal, [], exControl, net);
    figure; [mae, maeerr] = plotTwoSignals(si, S2);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);

    % simulate mVAR network with 1st frame & exogenous input signal
    netMVAR = initMvarNetwork(si, exSignal, [], exControl, lags);    
    [S3, time] = simulateMvarNetwork(si, exSignal, [], exControl, netMVAR);
    figure; [mae, maeerr] = plotTwoSignals(si, S3);

    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    figure; FC = plotFunctionalConnectivity(S3);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);
    figure; gcI = plotPairwiseGCI(S3);

    %% test pattern 2 
    % do training or load mVar DNN network
    lags = 5;
    mvardnnFile = ['results/mvardnn-sim-test8-4-' num2str(lags) '.mat'];
    if exist(mvardnnFile, 'file')
        load(mvardnnFile);
    else
        % init mVAR DNN network
        net = initDlcmNetwork(si, exSignal, [], exControl, lags);
        % training mVAR DNN network
        net = trainDlcmNetwork(si, exSignal, [], exControl, net, options);
        [time, loss, rsme] = getDlcmTrainingResult(net);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

        % init mVAR DNN network (linear)
        net2 = initDlcmNetwork(si, exSignal, [], exControl, lags, []);
        % training mVAR DNN network (linear)
        net2 = trainDlcmNetwork(si, exSignal, [], exControl, net2, options);
        [time, loss, rsme] = getDlcmTrainingResult(net2);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);

        save(mvardnnFile, 'net', 'net2');
    end

    % simulate mVAR DNN (DLCM) network with 1st frame & exogenous input signal
    [S, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, net);
    figure; [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);

    % simulate mVAR DNN (DLCM) linear network with 1st frame & exogenous input signal
    [S2, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, net2);
    figure; [mae, maeerr] = plotTwoSignals(si, S2);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);

    % simulate mVAR network with 1st frame & exogenous input signal
    netMVAR = initMvarNetwork(si, exSignal, [], exControl, lags);    
    [S3, time] = simulateMvarNetwork(si, exSignal, [], exControl, netMVAR);
    figure; [mae, maeerr] = plotTwoSignals(si, S3);

    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    figure; FC = plotFunctionalConnectivity(S2);
    figure; FC = plotFunctionalConnectivity(S3);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);
    figure; gcI = plotPairwiseGCI(S2);
    figure; gcI = plotPairwiseGCI(S3);
end

