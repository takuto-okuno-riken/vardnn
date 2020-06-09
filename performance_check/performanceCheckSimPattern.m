
function performanceCheckSimPattern
    % load signals
    load('test/testTrain-rand500-rand1.mat');
    siOrg = si;
    nodeNum = 8;
    inputNum = 4;
    sigLen = 100;
    si = siOrg(1:nodeNum,1:sigLen);
    inSignal = siOrg(nodeNum+1:nodeNum+inputNum,1:sigLen);
    % control is all positive input
    inControl = logical(ones(nodeNum,inputNum));

    %% pattern 1 -------------------------------------------------
%{
    disp('full random -- full independent nodes');
    si = siOrg(1:nodeNum,1:sigLen);
    inSignal = siOrg(nodeNum+1:nodeNum+inputNum,1:sigLen);
    inControl = logical(ones(nodeNum,inputNum));
    checkingPattern(si, inSignal, inControl, 1);
%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('node 2 and 6 are syncronized');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,:) = si(6,:);
    checkingPattern(si, inSignal, inControl, 2);
%%}
    %% pattern 3 -------------------------------------------------
%{
    disp('node 2 is excited by node 6');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = si(6,1:sigLen-1);
    checkingPattern(si, inSignal, inControl, 3);
%}
    %% pattern 4 -------------------------------------------------
%{
    disp('node 2 is excited half by node 6');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = si(6,1:sigLen-1) * 0.5;
    checkingPattern(si, inSignal, inControl, 4);
%}
    %% pattern 5 -------------------------------------------------
%%{
    disp('node 2,4 is excited by node 6');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);
    checkingPattern(si, inSignal, inControl, 5);
%%}
    %% pattern 6 -------------------------------------------------
%{
    disp('nodes are excited 6-.->2, 2-.->4');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,3:end) = si(6,1:sigLen-2);
    si(4,3:end) = si(2,1:sigLen-2);
    checkingPattern(si, inSignal, inControl, 6);
%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('node 2,4 is excited by exogenous input 2');
    si = siOrg(1:nodeNum, 1:sigLen);
    
    si(2,2:end) = inSignal(2,1:sigLen-1);
    si(4,2:end) = inSignal(2,1:sigLen-1);
    checkingPattern(si, inSignal, inControl, 7);
%%}
    %% pattern 8 -------------------------------------------------
%{
    disp('nodes are excited exIn2-.->2, 2-.->4');
    si = siOrg(1:nodeNum, 1:sigLen);
    si(2,3:end) = inSignal(2,1:sigLen-2);
    si(4,3:end) = si(2,1:sigLen-2);
    checkingPattern(si, inSignal, inControl, 8);
%}
end

function checkingPattern(si, inSignal, inControl, idx)
    nodeNum = size(si,1);
    inputNum = size(inSignal,1);
    sigLen = size(si,2);

    % do training or load DLCM network
    dlcmFile = ['performance_check/net-sim-pat' num2str(idx) '_' num2str(nodeNum) '-' num2str(inputNum) 'x' num2str(sigLen) '.mat'];
    if exist(dlcmFile, 'file')
        load(dlcmFile);
    else
        % init DLCM network
        netDLCM = initDlcmNetwork(si, inSignal, inControl);

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

        % training DLCM network
        netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
        [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
        disp(['train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
        %plotDlcmWeight(netDLCM);
        save(dlcmFile, 'netDLCM');
    end
    
    % simulate DLCM network with 1st frame & exogenous input signal
    [S, time] = simulateDlcmNetwork(si, inSignal, inControl, netDLCM);

    [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity([si; inSignal], inputNum);
    figure; FC = plotFunctionalConnectivity([S; inSignal], inputNum);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI([si; inSignal],3,10, inputNum);
    figure; gcI = plotPairwiseGCI([S; inSignal],3,10, inputNum);
    % show original time shifted correlation (tsc-FC)
    %tscFC = plotTimeShiftedCorrelation(si);
    % show deep-learning effective connectivity
    figure; dlEC = plotDlcmECmeanWeight(netDLCM);
end

