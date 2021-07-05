
function performanceCheckSimPattern
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

    %% pattern 1 -------------------------------------------------
%{
    disp('full random -- full independent nodes');
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    exControl = logical(ones(nodeNum,exNum));
    checkingPattern(si, exSignal, exControl, 1);
%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('node 2 and 6 are syncronized');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,:) = si(6,:);
    checkingPattern(si, exSignal, exControl, 2);
%%}
    %% pattern 3 -------------------------------------------------
%{
    disp('node 2 is excited by node 6');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = si(6,1:sigLen-1);
    checkingPattern(si, exSignal, exControl, 3);
%}
    %% pattern 4 -------------------------------------------------
%{
    disp('node 2 is excited half by node 6');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = si(6,1:sigLen-1) * 0.5;
    checkingPattern(si, exSignal, exControl, 4);
%}
    %% pattern 5 -------------------------------------------------
%%{
    disp('node 2,4 is excited by node 6');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);
    checkingPattern(si, exSignal, exControl, 5);
%%}
    %% pattern 6 -------------------------------------------------
%{
    disp('nodes are excited 6-.->2, 2-.->4');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,3:end) = si(6,1:sigLen-2);
    si(4,3:end) = si(2,1:sigLen-2);
    checkingPattern(si, exSignal, exControl, 6);
%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('node 2,4 is excited by exogenous input 2');
    si = siOrg(1:nodeNum, 1:sigLen);
    
    si(2,2:end) = exSignal(2,1:sigLen-1);
    si(4,2:end) = exSignal(2,1:sigLen-1);
    checkingPattern(si, exSignal, exControl, 7);
%%}
    %% pattern 8 -------------------------------------------------
%{
    disp('nodes are excited exIn2-.->2, 2-.->4');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,3:end) = exSignal(2,1:sigLen-2);
    si(4,3:end) = si(2,1:sigLen-2);
    checkingPattern(si, exSignal, exControl, 8);
%}
end

function checkingPattern(si, exSignal, exControl, idx)
    nodeNum = size(si,1);
    exNum = size(exSignal,1);
    sigLen = size(si,2);

    % do training or load VARDNN network
    netFile = ['results/net-sim-pat' num2str(idx) '_' num2str(nodeNum) '-' num2str(exNum) 'x' num2str(sigLen) '.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARDNN network
        netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);

        % set training options
        maxEpochs = 1000;
        sigLen = size(si,2);
        miniBatchSize = ceil(sigLen / 3);

        options = trainingOptions('adam', ...
            'ExecutionEnvironment','cpu', ...
            'MaxEpochs',maxEpochs, ...
            'MiniBatchSize',miniBatchSize, ...
            'Shuffle','every-epoch', ...
            'GradientThreshold',1,...
            'Verbose',false);
    %            'Plots','training-progress');

        % training VARDNN network
        netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
        [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        %plotMvarDnnWeight(netDLCM);
        save(netFile, 'netDLCM');
    end
    
    % simulate VARDNN network with 1st frame & exogenous input signal
    [S, time] = simulateMvarDnnNetwork(si, exSignal, [], exControl, netDLCM);

    [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si, exSignal, [], exControl, 1);
    figure; FC = plotFunctionalConnectivity(S, exSignal, [], exControl, 1);
    % show original & simulated signal granger causality index (GCI)
    figure; gcI = plotPairwiseGCI(si,exSignal,[],exControl,3,10,0.05,1);
    figure; gcI = plotPairwiseGCI(S,exSignal,[],exControl,3,10,0.05,1);
    % show original time shifted correlation (tsc-FC)
    %tscFC = plotTimeShiftedCorrelation(si);
    % show deep-learning effective connectivity
    figure; dlEC = plotMvarDnnDI(netDLCM,[],exControl,0,1);
end

