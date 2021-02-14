%%
% Return trained multivariate VAR DNN network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           mVAR DNN network structure
%  options       training options

function trainedNet = trainMvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options)
    global dlcmInitWeights;
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    trainedNet = net;
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end

    % training whole multivariate VAR DNN network
    disp('start training whole multivariate VAR DNN (DLCM) network');
    ticH = tic;
    nodeInputOrg = [];
    for i=1:lags, nodeInputOrg = [nodeInputOrg; X(:,i:end-(lags-i+1))]; end
    for i=1:lags, nodeInputOrg = [nodeInputOrg; exSignal(:,i:end-(lags-i+1))]; end
    
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,1);
    trainInfo = cell(nodeNum,1);
    initWeights = cell(nodeNum,1);
    for i=1:nodeNum
%    parfor i=1:nodeNum    % for parallel processing
        disp(['training node ' num2str(i)]);
        nodeTeach = X(i,1+lags:end);
        nodeInput = nodeInputOrg;
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', lags, size(nodeInput,2));
            nodeInput(1:nodeNum*lags,:) = nodeInput(1:nodeNum*lags,:) .* filter;
        end
        if ~isempty(exControl)
            filter = repmat(exControl(i,:).', lags, size(nodeInput,2));
            nodeInput(nodeNum*lags+1:end,:) = nodeInput(nodeNum*lags+1:end,:) .* filter;
        end
        [nodeNetwork{i}, trainInfo{i}] = trainNetwork(nodeInput, nodeTeach, nodeLayers{i}, options);
        initWeights{i} = dlcmInitWeights;
    end
    time = toc(ticH);
    trainedNet.nodeNum = nodeNum; % for compatibility
    trainedNet.exNum = exNum; % for compatibility
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.initWeights = initWeights;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole multivariate VAR DNN (DLCM) network! t = ' num2str(time) 's']);
end