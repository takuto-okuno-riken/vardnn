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
    Y = flipud(Y.'); % need to flip signal
    seqLen = sigLen-lags;
    
    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

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
        YTrain = Y(1:seqLen,i);
        iNodeControl = squeeze(nodeControl(i,:,:));
        iExControl = squeeze(exControl(i,:,:));
        if size(nodeControl,3)<=1, iNodeControl = iNodeControl.'; end
        if size(exControl,3)<=1, iExControl = iExControl.'; end
        for j=1:seqLen
            Xj = Y(j+1:j+lags,:).';
            Xj(1:nodeNum,:) = Xj(1:nodeNum,:) .* iNodeControl;
            Xj(nodeNum+1:end,:) = Xj(nodeNum+1:end,:) .* iExControl;
            XTrain{j} = Xj;
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