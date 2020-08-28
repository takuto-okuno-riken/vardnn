
function testSimulation
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
    dlcmFile = ['results/dlcm-gc-test' num2str(nodeNum) '-' num2str(inputNum) '.mat'];
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
        [netDLCM, time] = recoveryTrainDlcmNetwork(si, inSignal, [], inControl, netDLCM, options);
        save(dlcmFile, 'netDLCM');
    end
    
    % simulate DLCM network with 1st frame & exogenous input signal
    [S, time] = simulateDlcmNetwork(si, inSignal, [], inControl, netDLCM);

    [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);

    % show DLCM-GC
    figure; dlGC = plotDlcmGCI(si, inSignal, [], inControl, netDLCM);
end

