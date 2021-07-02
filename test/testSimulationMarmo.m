
function testSimulationMarmo
    % load signals
    load('data/marmoset-aneth-sample1-roi52.mat');
%    figure; FC = plotFunctionalConnectivity(si);
%    figure; gcI = plotPairwiseGCI(si);

    siOrg = si;
    nodeNum = 12;
    exNum = 10;
    sigLen = 300;
    si = bold2dnnSignal(siOrg(1:nodeNum,1:sigLen), 0.2);
    exSignal = bold2dnnSignal(siOrg(nodeNum+1:nodeNum+exNum,1:sigLen), 0.2);
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
        'Shuffle','every-epoch', ...
        'GradientThreshold',5,...
        'L2Regularization',0.1, ... % for gaussian distribution (Marmo)
        'Verbose',false);
%            'Plots','training-progress');

    %% test pattern 1 
    % do training or load DLCM network
    dlcmFile = ['results/dlcm-sim-marmo' num2str(nodeNum) '-' num2str(exNum) '.mat'];
    if exist(dlcmFile, 'file')
        load(dlcmFile);
    else
        % init VARDNN network
        netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
        % training VARDNN network
        netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
        % recover training 
        [netDLCM, time] = recoveryTrainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
        [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);        save(dlcmFile, 'netDLCM');
    end
    
    % simulate DLCM network with 1st frame & exogenous input signal
    [S, time] = simulateDlcmNetwork(si, exSignal, [], exControl, netDLCM);

    figure; [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
%    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si,[],[],[],3,0);
%    figure; gcI = plotPairwiseGCI(S,[],[],[],3,0);
    % show DLCM-GC
    figure; dlGC = plotMvarDnnGCI(si, exSignal, [], exControl, netDLCM, 0);
end

