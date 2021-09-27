%%
% Return trained multivariate PC VAR DNN network
% input :
%  X             multivariate time series matrix (node x time series)
%  pcX           multivariate time series matrix (pcNode x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           mVAR DNN network structure
%  options       training options

function trainedNet = trainMpcvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options)
    global dlcmInitWeights;
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    trainedNet = net;
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    inputNum = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal
    
    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % set input signals
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end

    % training whole multivariate PC VAR DNN network
    disp('start training whole multivariate PC VAR DNN network');
    ticH = tic;
    
    mu = net.mu;
    coeff = net.coeff;
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,1);
    trainInfo = cell(nodeNum,1);
    initWeights = cell(nodeNum,1);

    % train DNN network
    for i=1:nodeNum
%    parfor i=1:nodeNum    % for parallel processing
        disp(['training node ' num2str(i)]);
        [~,idx] = find(control(i,:,:)==1);
        
        % VARDNN teaching signals
        Xt = Y(1:sigLen-lags,i).';
        Xti = Yj(:,idx);
        score = (Xti - mu{i}) / coeff{i}.'; % calculate score
        
        [nodeNetwork{i}, trainInfo{i}] = trainNetwork(score.', Xt, nodeLayers{i}, options);
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
    disp(['finish training whole multivariate PC VAR DNN network! t = ' num2str(time) 's']);
end