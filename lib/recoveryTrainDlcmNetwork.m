%%
% Return trained DLCM network
% input :
%  X             multivariate time series matrix (node x time series)
%  inSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  inControl     exogenous input control matrix for each node (node x exogenous input)
%  net           DLCM network structure
%  options       training options
%  maeTh         MAE threshold to finish recovery training (default:0.05)

function [trainedNet, time, mae] = recoveryTrainDlcmNetwork(X, inSignal, nodeControl, inControl, net, options, maeTh)
    if nargin < 7, maeTh = 0.05; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    trainedNet = net;

    ticH = tic; % start stop watch

    % simulate DLCM network with 1st frame & exogenous input signal
    [S, time] = simulateDlcmNetwork(X, inSignal, nodeControl, inControl, net);
    [mae, ~] = getTwoSignalsError(X, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    if mae < maeTh
        time = toc(ticH); return;
    end

    % recovery training whole DLCM network
    lastsi = X(:,1:end-1);
    lastinSignal = inSignal(:,1:end-1);
    lastTeach = X(:,2:end);
    lastS = S;
    lastmae = mae;
%    correct = 0;
    rtes = [10, 20, 30, ...
        15, 15+floor(sigLen*0.1), 15+floor(sigLen*0.2), ...
        20, 20+floor(sigLen*0.15), 20+floor(sigLen*0.3)];
    maxEpochs = [400, 600, 800, 1000, 1200, 1200, 800, 1000, 1000];

% for debug ---------------
%{
dlcmFile = ['results/net-sim-err-inRecov_' num2str(nodeNum) 'x' num2str(sigLen) '.mat'];
if exist(dlcmFile, 'file')
    load(dlcmFile);
else
    nets = {};
end
%}
% for debug ---------------
    disp('start recovery training whole DLCM network');
    for k=1:length(rtes)
%{
        rts = 6 + (k-1)*10;
        rte = 20 + (k-1)*10;
        if isempty(inSignal)
            nodeInput = [lastsi, S(:,rts:rte)];
            lastsi = nodeInput;
        else
            si2 = [lastsi, S(:,rts:rte)];
            inSig2 = [lastinSignal, inSignal(:,rts:rte)];
            nodeInput = [si2; inSig2];
            lastsi = si2;
            lastinSignal = inSig2;
        end
        nodeTeach = [lastTeach, si(:,rts+1:rte+1)]; 
        lastTeach = nodeTeach;
%}
        rte = rtes(k);
%        rte = 10 + floor((k-1)*(sigLen/10) + correct);
        if isempty(inSignal)
            nodeInputOrg = [lastsi, lastS(:,1:rte)];
        else
            si2 = [lastsi, lastS(:,1:rte)];
            inSig2 = [lastinSignal, inSignal(:,1:rte)];
            nodeInputOrg = [si2; inSig2];
        end
        nodeTeach = [lastTeach, X(:,2:rte+1)]; 
%{
        rte = 20 + (k-1)*10;
%        rte = 99;
        if isempty(inSignal)
            nodeInput = S(:,1:rte);
        else
            si2 = S(:,1:rte);
            inSig2 = inSignal(:,1:rte);
            nodeInput = [si2; inSig2];
        end
        nodeTeach = si(:,2:rte+1); 
%}
        % learning option
        trainOpt = trainingOptions('adam', ...
            'ExecutionEnvironment',options.ExecutionEnvironment, ...
            'MaxEpochs',maxEpochs(k), ...
            'MiniBatchSize',size(nodeTeach,2), ...
            'Shuffle',options.Shuffle, ...
            'L2Regularization',options.L2Regularization, ...
            'GradientThreshold',options.GradientThreshold, ...
            'GradientThresholdMethod', options.GradientThresholdMethod, ...
            'Verbose',options.Verbose, ...
            'Plots',options.Plots);
% for debug ---------------
%{
if length(nets) >= k
    net = nets{k}; % for debug
else
%}
% for debug ---------------
        ticH2 = tic;
        for i=1:nodeNum
            disp(['training node ' num2str(i)]);
            nodeInput = nodeInputOrg;
            if ~isempty(nodeControl)
                filter = repmat(nodeControl(i,:).', 1, size(nodeInput,2));
                nodeInput(1:nodeNum,1) = nodeInput(1:nodeNum,1) .* filter;
            end
            if ~isempty(inControl)
                filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
                nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
            end
            [net.nodeNetwork{i}, net.trainInfo{i}] = trainNetwork(nodeInput, nodeTeach(i,:), trainedNet.nodeNetwork{i}.Layers, trainOpt);
        end
        time = toc(ticH2);
        disp(['recovery training (' num2str(k) '): train time=' num2str(time)]);
% for debug ---------------
%{
nets{end+1} = net; % for debug
save(dlcmFile, 'nets');
end
%}
% for debug ---------------
        % simulate again DLCM network with 1st frame & exogenous input signal
        [S, time] = simulateDlcmNetwork(X, inSignal, nodeControl, inControl, net);
        [err, ~] = getTwoSignalsError(X, S);
        disp(['recovery training (' num2str(k) '): simulation time=' num2str(time) ', mae=' num2str(err)]);
%{
if exist('f')>0, figure(f); else f = figure;  end
plotTwoSignals(X, S);
drawnow;
%}
        if lastmae > err
            lastmae = err; mae = err;
            lastS = S;
            trainedNet.nodeNetwork = net.nodeNetwork;
            trainedNet.trainInfo = net.trainInfo;
            if mae < maeTh
                k = length(rtes); % finish loop
            end
        else
%            correct = floor(correct - (sigLen/10) * 3 / 4);
        end
    end
    time = toc(ticH);
    trainedNet.recoverTrainTime = time;
    disp(['finish recovery training whole DLCM network! t = ' num2str(time) 's']);
end