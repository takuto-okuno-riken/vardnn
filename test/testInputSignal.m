
function testInputSignal
    % load signals
    load('test/testTrain-rand500-rand1.mat');
    siOrg = si;
    nodeNum = 8;
    inputNum = 4;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    inSignal = siOrg(nodeNum+1:nodeNum+inputNum,1:sigLen);

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

    %% test pattern 1 -- no input signal (default weight initializer)
%%{
    % init DLCM network
    netDLCM = initDlcmNetwork(si, []);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, [], [], netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['1) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
%%}
    %% test pattern 2 -- input signal without input control (default weight initializer)
%{
    % init DLCM network
    netDLCM = initDlcmNetwork(si, inSignal, []);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, inSignal, [], netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['2) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%}
    %% test pattern 3 -- input signal with input control (default weight initializer)
%%{
    % control is all zero
    inControl = logical(zeros(nodeNum,inputNum));
    % init DLCM network
    netDLCM = initDlcmNetwork(si, inSignal, inControl);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%%}
    %% test pattern 4 -- input signal with input control (default weight initializer)
%{
    % control one by one
    inControl = logical(zeros(nodeNum,inputNum));
    for i=1:nodeNum
        inControl(i, mod(i-1,4)+1) = 1;
    end
    % init DLCM network
    netDLCM = initDlcmNetwork(si, inSignal, inControl);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%}
    %% test pattern 5 -- input signal with input control (default weight initializer)
%%{
    % control one by one
    inControl = logical(zeros(nodeNum,inputNum));
    for i=1:nodeNum
        inControl(i, mod(i-1,4)+1) = 1;
        inControl(i, mod(i+1,4)+1) = 1;
    end
    % init DLCM network
    netDLCM = initDlcmNetwork(si, inSignal, inControl);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%%}
    %% test pattern 6 -- input signal with input control (default weight initializer)
%{
    % control one by one
    inControl = logical(ones(nodeNum,inputNum));
    % init DLCM network
    netDLCM = initDlcmNetwork(si, inSignal, inControl);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%}
    %% test pattern 7 -- input signal with input control without weight initializer
%%{
    % control one by one
    inControl = logical(zeros(nodeNum,inputNum));
    for i=1:nodeNum
        inControl(i, mod(i-1,4)+1) = 1;
        inControl(i, mod(i+1,4)+1) = 1;
    end
    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(nodeNum, sigLen);
    % layer parameters
    netDLCM = createDlcmNetwork(nodeNum, inputNum, hiddenNums, inControl, []);
    % training DLCM network
    netDLCM = trainDlcmNetwork(si, inSignal, inControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%%}
end

