%%
% Return trained DLCM network
% input :
%  X             multivariate time series matrix (node x time series)
%  inSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  inControl     exogenous input control matrix for each node (node x exogenous input)
%  net           DLCM network structure
%  options       training options

function trainedNet = trainDlcmNetwork(X, inSignal, inControl, net, options)
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
    for i=1:nodeNum
        disp(['training node ' num2str(i)]);
        nodeTeach = X(i,2:end);
        nodeInput = nodeInputOrg;
        if ~isempty(inControl)
            filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        [trainedNet.nodeNetwork{i}, trainedNet.trainInfo{i}] = trainNetwork(nodeInput, nodeTeach, trainedNet.nodeLayers{i}, options);
        trainedNet.initWeights{i} = dlcmInitWeights;
    end
    time = toc(ticH);
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole DLCM network! t = ' num2str(time) 's']);
end