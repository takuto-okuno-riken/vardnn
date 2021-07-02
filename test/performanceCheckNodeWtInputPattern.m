
function performanceCheckNodeWtInputPattern
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    
    nodeNum = 8;
    exNum = 2;
    sigLen = 100;

    si = siOrg(1:nodeNum, 1:sigLen);
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    exControl = logical(ones(nodeNum,exNum));

    %% pattern 1 -------------------------------------------------
%{
    disp('full random -- full independent nodes');
    checkingPattern(si, exSignal, exControl, 1);
%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('node 2 and exogenous input1 are syncronized');
    si(2,:) = exSignal(1,:);
    checkingPattern(si, exSignal, exControl, 2);
%%}
    %% pattern 3 -------------------------------------------------
%%{
    disp('node 2 is excited by exogenous input 1');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = exSignal(1,1:sigLen-1);
    checkingPattern(si, exSignal, exControl, 3);
%%}
    %% pattern 4 -------------------------------------------------
%{
    disp('node 2 is excited half by exogenous input 1');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = exSignal(1,1:sigLen-1) * 0.5;
    checkingPattern(si, exSignal, exControl, 4);
%}
    %% pattern 5 -------------------------------------------------
%%{
    disp('node 2,4 is excited by exogenous input 2');
    si = siOrg(1:nodeNum, 1:sigLen);
    
    si(2,2:end) = exSignal(2,1:sigLen-1);
    si(4,2:end) = exSignal(2,1:sigLen-1);
    checkingPattern(si, exSignal, exControl, 5);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('nodes are excited exIn2-.->2, 2-.->4');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = exSignal(2,1:sigLen-1);
    si(4,2:end) = si(2,1:sigLen-1);
    checkingPattern(si, exSignal, exControl, 6);
%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('node 2,4 and exogenous input1 are syncronized, only 4 receive input1');
    si = siOrg(1:nodeNum, 1:sigLen);
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    si(2,:) = exSignal(1,:);
    si(4,:) = exSignal(1,:);
    exControl = logical(zeros(nodeNum,exNum));
    exControl(4,1) = 1;
    checkingPattern(si, exSignal, exControl, 7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('node 2,4 are excited by exogenous input1, only 4 receive input1');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = exSignal(1,1:sigLen-1);
    si(4,3:end) = exSignal(1,2:sigLen-1);
    exControl = logical(zeros(nodeNum,exNum));
    exControl(4,1) = 1;
    checkingPattern(si, exSignal, exControl, 8);
%%}
end

%% 
function [FC, dlEC, gcI] = checkingPattern(si, exSignal, exControl, idx)
    nodeNum = size(si,1);
    exNum = size(exSignal,1);
    sigLen = size(si,2);

    % layer parameters
    netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);

    % show signals before training
    %{
    maxEpochs = 1;
    miniBatchSize = 1;
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'Shuffle','every-epoch', ...
        'GradientThreshold',1,...
        'Verbose',false);
%            'Plots','training-progress');

    disp('initial state before training');
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    [t,mae,maeerr] = plotNodeSignals(nodeNum,si,exSignal,netDLCM);
    disp(['t=' num2str(t) ', mae=' num2str(mae)]);
    %}
    % training VARDNN network
    maxEpochs = 1000;
    miniBatchSize = ceil(sigLen / 3);
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'Shuffle','every-epoch', ...
        'GradientThreshold',1,...
        'Verbose',false);
%            'Plots','training-progress');

    disp('start training');
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    dlcmFile = ['results/net-pat-in' num2str(idx) '.mat'];
    save(dlcmFile, 'netDLCM');

    % show signals after training
    [S, t,mae,maeerr] = plotPredictSignals(si,exSignal,[],exControl,netDLCM);
    disp(['t=' num2str(t) ', mae=' num2str(mae)]);

    % show original signal FC
    figure; FC = plotFunctionalConnectivity(si,exSignal,[],exControl,1);
    % show original signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si,exSignal,[],exControl,3,10,0.05,1);
    % show original time shifted correlation (tsc-FC)
    %figure; tscFC = plotTimeShiftedCorrelation(si);
    % show deep-learning effective connectivity
    figure; dlEC = plotMvarDnnEC(netDLCM, [], exControl, 0, 1);    
end

