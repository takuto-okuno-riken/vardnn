
function testInputSignal
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 4;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);

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

    %% test pattern 1 -- no exogenous signal (default weight initializer)
%%{
    % init VARDNN network
    netDLCM = initMvarDnnNetwork(si);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, [], [], [], netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['1) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
%%}
    %% test pattern 2 -- exogenous signal without exogenous control (default weight initializer)
%{
    % init VARDNN network
    netDLCM = initMvarDnnNetwork(si, exSignal, []);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], [], netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['2) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%}
    %% test pattern 3 -- exogenous signal with exogenous control (default weight initializer)
%%{
    % control is all zero
    exControl = logical(zeros(nodeNum,exNum));
    % init VARDNN network
    netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%%}
    %% test pattern 4 -- exogenous signal with exogenous control (default weight initializer)
%{
    % control one by one
    exControl = logical(zeros(nodeNum,exNum));
    for i=1:nodeNum
        exControl(i, mod(i-1,4)+1) = 1;
    end
    % init VARDNN network
    netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%}
    %% test pattern 5 -- exogenous signal with exogenous control (default weight initializer)
%%{
    % control one by one
    exControl = logical(zeros(nodeNum,exNum));
    for i=1:nodeNum
        exControl(i, mod(i-1,4)+1) = 1;
        exControl(i, mod(i+1,4)+1) = 1;
    end
    % init VARDNN network
    netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%%}
    %% test pattern 6 -- exogenous signal with exogenous control (default weight initializer)
%{
    % control one by one
    exControl = logical(ones(nodeNum,exNum));
    % init VARDNN network
    netDLCM = initMvarDnnNetwork(si, exSignal, [], exControl);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%}
    %% test pattern 7 -- exogenous signal with exogenous control without weight initializer
%%{
    % control one by one
    exControl = logical(zeros(nodeNum,exNum));
    for i=1:nodeNum
        exControl(i, mod(i-1,4)+1) = 1;
        exControl(i, mod(i+1,4)+1) = 1;
    end
    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(nodeNum, sigLen);
    % layer parameters
    netDLCM = createMvarDnnNetwork(nodeNum, exNum, hiddenNums, 1, [], exControl, []);
    % training VARDNN network
    netDLCM = trainMvarDnnNetwork(si, exSignal, [], exControl, netDLCM, options);
    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
    disp(['3) train result time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    plotDlcmWeight(netDLCM);
%%}
end

