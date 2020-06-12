
function performanceCheckNodeWtInputPattern
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    
    nodeNum = 8;
    inputNum = 2;
    sigLen = 100;

    si = siOrg(1:nodeNum, 1:sigLen);
    inSignal = siOrg(nodeNum+1:nodeNum+inputNum,1:sigLen);
    inControl = logical(ones(nodeNum,inputNum));
    
    %% pattern 1 -------------------------------------------------
%{
    disp('full random -- full independent nodes');
    checkingPattern(si, inSignal, inControl, 1);
%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('node 2 and exogenous input1 are syncronized');
    si(2,:) = inSignal(1,:);
    checkingPattern(si, inSignal, inControl, 2);
%%}
    %% pattern 3 -------------------------------------------------
%%{
    disp('node 2 is excited by exogenous input 1');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = inSignal(1,1:sigLen-1);
    checkingPattern(si, inSignal, inControl, 3);
%%}
    %% pattern 4 -------------------------------------------------
%{
    disp('node 2 is excited half by exogenous input 1');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = inSignal(1,1:sigLen-1) * 0.5;
    checkingPattern(si, inSignal, inControl, 4);
%}
    %% pattern 5 -------------------------------------------------
%%{
    disp('node 2,4 is excited by exogenous input 2');
    si = siOrg(1:nodeNum, 1:sigLen);
    
    si(2,2:end) = inSignal(2,1:sigLen-1);
    si(4,2:end) = inSignal(2,1:sigLen-1);
    checkingPattern(si, inSignal, inControl, 5);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('nodes are excited exIn2-.->2, 2-.->4');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = inSignal(2,1:sigLen-1);
    si(4,2:end) = si(2,1:sigLen-1);
    checkingPattern(si, inSignal, inControl, 6);
%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('node 2,4 and exogenous input1 are syncronized, only 4 receive input1');
    si = siOrg(1:nodeNum, 1:sigLen);
    inSignal = siOrg(nodeNum+1:nodeNum+inputNum,1:sigLen);
    si(2,:) = inSignal(1,:);
    si(4,:) = inSignal(1,:);
    inControl = logical(zeros(nodeNum,inputNum));
    inControl(4,1) = 1;
    checkingPattern(si, inSignal, inControl, 7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('node 2,4 are excited by exogenous input1, only 4 receive input1');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = inSignal(1,1:sigLen-1);
    si(4,3:end) = inSignal(1,2:sigLen-1);
    inControl = logical(zeros(nodeNum,inputNum));
    inControl(4,1) = 1;
    checkingPattern(si, inSignal, inControl, 8);
%%}
end

%% 
function [FC, dlEC, gcI] = checkingPattern(si, inSignal, inControl, idx)
    nodeNum = size(si,1);
    inputNum = size(inSignal,1);
    sigLen = size(si,2);

    % layer parameters
    netDLCM = initDlcmNetwork(si, inSignal, inControl);

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
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    [t,mae,maeerr] = plotNodeSignals(nodeNum,si,inSignal,netDLCM);
    disp(['t=' num2str(t) ', mae=' num2str(mae)]);
    %}
    % training DLCM network
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
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    dlcmFile = ['performance_check/net-pat-' num2str(idx) '.mat'];
    save(dlcmFile, 'netDLCM');

    % show signals after training
    [S, t,mae,maeerr] = plotPredictSignals(si,inSignal,inControl,netDLCM);
    disp(['t=' num2str(t) ', mae=' num2str(mae)]);

    siWithInput = [si; inSignal];
    % show original signal FC
    FC = plotFunctionalConnectivity(siWithInput,inputNum);
    % show original signal granger causality index (gc-EC)
    gcI = plotPairwiseGCI(siWithInput,3,10,inputNum);
    % show original time shifted correlation (tsc-FC)
    %tscFC = plotTimeShiftedCorrelation(si);
    % show deep-learning effective connectivity
    dlEC = plotDlcmECmeanWeight(netDLCM);    
end

