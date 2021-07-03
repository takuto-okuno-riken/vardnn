%%
% Return trained multivariate VAR DNN network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input)
%  net           mVAR DNN network structure
%  options       training options
%  maeMin        minimum MAE of re-training (default:0.008)
%  retrainRange  error range rate(0-1) of re-training (default:[0.1, 0.3])
%  retrainMax    muximum re-training number (default:10)
%  retrainEpochs DNN retraining epoch number (default:400)

function [trainedNet, time, mae, errs] = recoveryTrainMvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options, maeMin, retrainRange, retrainMax, retrainEpochs)
    if nargin < 10, retrainEpochs = 400; end
    if nargin < 9, retrainMax = 10; end
    if nargin < 8, retrainRange = [0.1, 0.3]; end
    if nargin < 7, maeMin = 0.008; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    trainedNet = net;

    ticH = tic; % start stop watch

    % simulate mVAR DNN network with 1st frame & exogenous input signal
    [S, time] = simulateMvarDnnNetwork(X, exSignal, nodeControl, exControl, net);
    [mae, ~, errs] = getTwoSignalsError(X, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);

    errX = mean(abs(errs),1);
    errMax = max(errX);
    errRange(1) = retrainRange(1) * errMax;
    errRange(2) = retrainRange(2) * errMax;
%    fg = figure; plot(errX);

    errIdx = find(errRange(1) < errX & errX < errRange(2));
    if mae < maeMin
        time = toc(ticH); return;
    end
    if isempty(errIdx) % bad case
        [B, I] = sort(errX);
        for i=1:length(B)
            if B(i)>0, errIdx = [errIdx, I(i)]; end
            if length(errIdx)>2, break; end
        end
    end
    lastmae = mae;

    % recovery training whole VARDNN network
    baseS = [];
    for i=1:lags, baseS = [baseS; X(:,i:end-(lags-i+1))]; end
    for i=1:lags, baseS = [baseS; exSignal(:,i:end-(lags-i+1))]; end
    baseTeach = X(:,lags+1:end);
    lastS = S;
    inputHist = [];
    teachHist = [];

    disp('start recovery training whole VARDNN network');
    for k=1:retrainMax
        if errIdx(end)==sigLen, errIdx(end)=[]; end
        errInput = [];
        for j=1:length(errIdx)
            eY = [];
            for i=1:lags, eY = [eY; lastS(:,errIdx(j)-(i-1))]; end
            if ~isempty(exSignal)
                for i=1:lags, eY = [eY; exSignal(:,errIdx(j)-(i-1))]; end
            end
            errInput = [errInput, eY];
        end
        errTeach = X(:,errIdx+1);
        nodeInputOrg = [baseS, inputHist, errInput];
        nodeTeach = [baseTeach, teachHist, errTeach];

        % learning option
        trainOpt = trainingOptions('adam', ...
            'ExecutionEnvironment',options.ExecutionEnvironment, ...
            'MaxEpochs',retrainEpochs, ...
            'MiniBatchSize',size(nodeTeach,2), ...
            'Shuffle',options.Shuffle, ...
            'L2Regularization',options.L2Regularization, ...
            'GradientThreshold',options.GradientThreshold, ...
            'GradientThresholdMethod', options.GradientThresholdMethod, ...
            'Verbose',options.Verbose, ...
            'Plots',options.Plots);

        ticH2 = tic;
        nodeNetwork = trainedNet.nodeNetwork;
        trainInfo = cell(nodeNum,1);
%        for i=1:nodeNum
        parfor i=1:nodeNum
            disp(['training node ' num2str(i)]);
            nodeInput = nodeInputOrg;
            if ~isempty(nodeControl)
                filter = repmat(nodeControl(i,:).', lags, size(nodeInput,2));
                nodeInput(1:nodeNum*lags,:) = nodeInput(1:nodeNum*lags,:) .* filter;
            end
            if ~isempty(exControl)
                filter = repmat(exControl(i,:).', lags, size(nodeInput,2));
                nodeInput(nodeNum*lags+1:end,:) = nodeInput(nodeNum*lags+1:end,:) .* filter;
            end
            [nodeNetwork{i}, trainInfo{i}] = trainNetwork(nodeInput, nodeTeach(i,:), nodeNetwork{i}.Layers, trainOpt);
        end
        net.nodeNetwork = nodeNetwork;
        net.trainInfo = trainInfo;
        time = toc(ticH2);
        disp(['recovery training (' num2str(k) '): train time=' num2str(time)]);

        % simulate again VARDNN network with 1st frame & exogenous input signal
        [S, time] = simulateMvarDnnNetwork(X, exSignal, nodeControl, exControl, net);
        [mae, ~, errs2] = getTwoSignalsError(X, S);
        disp(['recovery training (' num2str(k) '): simulation time=' num2str(time) ', mae=' num2str(mae)]);

        errX = mean(abs(errs2),1);
        errs = cat(3, errs, errs2);
        errIdx = find(errRange(1) < errX & errX < errRange(2));
        if lastmae > mae
            lastmae = mae;
            lastS = S;
            inputHist = [inputHist, errInput];
            teachHist = [teachHist, errTeach];
            trainedNet.nodeNetwork = nodeNetwork;
            trainedNet.trainInfo = trainInfo;
        end
        if mae < maeMin
            break;
        end
        if isempty(errIdx) % bad case
            [B, I] = sort(errX);
            for i=1:length(B)
                if B(i)>0, errIdx = [errIdx, I(i)]; end
                if length(errIdx)>2, break; end
            end
        end
    end
    time = toc(ticH);
    trainedNet.recoverTrainTime = time;
    disp(['finish recovery training whole VARDNN network! t = ' num2str(time) 's']);
end