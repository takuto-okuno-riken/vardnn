%%
% Return trained DLCM network
% input :
%  X             multivariate time series matrix (node x time series)
%  inSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  inControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           DLCM network structure
%  options       training options

function trainedNet = trainDlcmNetwork(X, inSignal, nodeControl, inControl, net, options)
    global dlcmInitWeights;
    nodeNum = length(net.nodeLayers);
    trainedNet = net;

    % training whole DLCM network
    disp('start training whole DLCM network');
    ticH = tic;
    if isempty(inSignal)
        nodeInputOrg = X(:,1:end-1);
    else
        nodeInputOrg = [X(:,1:end-1); inSignal(:,1:end-1)];
    end
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,1);
    trainInfo = cell(nodeNum,1);
    initWeights = cell(nodeNum,1);
    for i=1:nodeNum
%    parfor i=1:nodeNum    % for parallel processing
        disp(['training node ' num2str(i)]);
        nodeTeach = X(i,2:end);
        nodeInput = nodeInputOrg;
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', 1, size(nodeInput,2));
            nodeInput(1:nodeNum,:) = nodeInput(1:nodeNum,:) .* filter;
        end
        if ~isempty(inControl)
            filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        [nodeNetwork{i}, trainInfo{i}] = trainNetwork(nodeInput, nodeTeach, nodeLayers{i}, options);
        initWeights{i} = dlcmInitWeights;
    end
    time = toc(ticH);
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.initWeights = initWeights;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole DLCM network! t = ' num2str(time) 's']);
end