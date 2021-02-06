%%
% Return trained pairwised VAR DNN network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           Pairwised VAR DNN network structure
%  options       training options

function trainedNet = trainPvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options)
    global dlcmInitWeights;
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    sigLen = size(X,2);
    trainedNet = net;
    lags = net.lags;
    
    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;

    % training whole Pairwised VAR DNN network
    disp('start training whole pairwised VAR DNN network');
    ticH = tic;
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,nodeMax);
    trainInfo = cell(nodeNum,nodeMax);
    initWeights = cell(nodeNum,nodeMax);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            disp(['training node ' num2str(i) '-' num2str(j)]);
            nodeTeach = Y(i,lags+1:end);
            nodeInput = nan(2*lags,sigLen-lags);
            for k=1:lags, nodeInput(k,:) = Y(i,k:end-lags+(k-1)); end
            for k=1:lags, nodeInput(lags+k,:) = Y(j,k:end-lags+(k-1)); end
            [nodeNetwork{i,j}, trainInfo{i,j}] = trainNetwork(nodeInput, nodeTeach, nodeLayers{i,j}, options);
            initWeights{i,j} = dlcmInitWeights;
        end
    end
    time = toc(ticH);
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.initWeights = initWeights;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole Pairwised VAR DNN network! t = ' num2str(time) 's']);
end