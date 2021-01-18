
function testRecoverTrain
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

    % set training options
    maxEpochs = 1000;
    sigLen = size(si,2);
    miniBatchSize = ceil(sigLen / 3);

    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'L2Regularization',0.01, ...
        'Shuffle','every-epoch', ...
        'GradientThreshold',1,...
        'Verbose',false);
%            'Plots','training-progress');

    %% test pattern 1 
    % do training or load DLCM network
    dlcmFile = 'results/dlcm-sim-test8-4.mat';
    if exist(dlcmFile, 'file')
        load(dlcmFile);
    else
        % init DLCM network
        netDLCM = initDlcmNetwork(si, exSignal, [], exControl);
        % training DLCM network
        netDLCM = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
        [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        %plotDlcmWeight(netDLCM);
        save(dlcmFile, 'netDLCM');
    end
    
    % recoverty training
    [netDLCM, time] = recoveryTrainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);

    % show result of recoverty training
    [S, time] = simulateDlcmNetwork(si, exSignal, [], exControl, netDLCM);
    [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    save(dlcmFile, 'netDLCM');
%%{
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);
    % show original time shifted correlation (tsc-FC)
    %tscFC = plotTimeShiftedCorrelation(si);
    % show deep-learning effective connectivity
    figure; dlEC = plotDlcmECmeanWeight(netDLCM);
    % plot correlation graph between original predicted node signals
    figure; R = plotTwoSignalsCorrelation(si, S) % show R result
    figure; R = getCosSimilarity(si, S) % show R result
%%}
end

