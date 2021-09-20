%%
% Return trained multivariate VAR LSTM network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           mVAR LSTM network structure
%  options       training options

function trainedNet = trainMvarLstmNetwork(X, exSignal, nodeControl, exControl, net, options)
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    sigLen = size(X,2);
    trainedNet = net;
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end

    % set node input
    Y = [X; exSignal];
    seqLen = sigLen-lags;

    % training whole multivariate VAR LSTM network
    disp('start training whole multivariate VAR LSTM network');
    ticH = tic;
    
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,1);
    trainInfo = cell(nodeNum,1);
    initWeights = cell(nodeNum,1);
    for i=1:nodeNum
%    parfor i=1:nodeNum    % for parallel processing
        disp(['training node ' num2str(i)]);
        XTrain = cell(seqLen,1);
        YTrain = X(i,1+lags:end).';
        nodeInput = Y;
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', 1, size(nodeInput,2));
            nodeInput(1:nodeNum,:) = nodeInput(1:nodeNum,:) .* filter;
        end
        if ~isempty(exControl)
            filter = repmat(exControl(i,:).', 1, size(nodeInput,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        for j=1:seqLen
            XTrain{j} = nodeInput(:,j:j+(lags-1));
        end
        [nodeNetwork{i}, trainInfo{i}] = trainNetwork(XTrain, YTrain, nodeLayers{i}, options);
    end
    time = toc(ticH);
    trainedNet.nodeNum = nodeNum; % for compatibility
    trainedNet.exNum = exNum; % for compatibility
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.initWeights = initWeights;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole multivariate VAR LSTM network! t = ' num2str(time) 's']);
end