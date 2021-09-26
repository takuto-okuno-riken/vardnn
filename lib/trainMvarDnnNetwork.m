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
    global dnnInitWeights;
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    trainedNet = net;
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end

    % set node input
    Y = [X; exSignal];
    inputNum = nodeNum + exNum;

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    Y = flipud(Y.'); % need to flip signal

    % training whole multivariate VAR DNN network
    disp('start training whole multivariate VAR DNN network');
    ticH = tic;

    % set training data
    Yj = zeros(sigLen-lags, lags*inputNum);
    for p=1:lags
        Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sigLen-lags+p,:);
    end
    
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,1);
    trainInfo = cell(nodeNum,1);
    initWeights = cell(nodeNum,1);
    for i=1:nodeNum
%    parfor i=1:nodeNum    % for parallel processing
        disp(['training node ' num2str(i)]);
        [~,idx] = find(control(i,:,:)==1);
        
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);

        % training network
        [nodeNetwork{i}, trainInfo{i}] = trainNetwork(Xti.', Xt.', nodeLayers{i}, options);
        initWeights{i} = dnnInitWeights;
    end
    time = toc(ticH);
    trainedNet.nodeNum = nodeNum; % for compatibility
    trainedNet.exNum = exNum; % for compatibility
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.initWeights = initWeights;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole multivariate VAR DNN network! t = ' num2str(time) 's']);
end