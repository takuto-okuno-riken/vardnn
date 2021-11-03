%%
% Return trained Nonlinear VAR Neural Network (NVARNN) network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (1 x node) (optional)
%  exControl     exogenous input control matrix for each node (1 x exogenous input) (optional)
%  net           nVARNN network structure
%  options       training options

function trainedNet = trainNvarnnNetwork(X, exSignal, nodeControl, exControl, net, options)
    global dnnInitWeights;
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    trainedNet = net;
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    inputNum = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal
    
    % set control 3D matrix (1 x node x lags)
    if isempty(nodeControl)
        nodeControl = ones(1,nodeNum,lags);
    end
    if isempty(exControl)
        exControl = ones(1,exNum,lags);
    end
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    cIdx = find(control(:)==1);

    % set input signals
    Yj = zeros(sigLen-lags, lags*inputNum);
    for p=1:lags
        Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sigLen-lags+p,:);
    end

    % training whole nVARNN network
    disp('start training whole nVARNN network');
    ticH = tic;

    Xt = Y(1:sigLen-lags,1:nodeNum);
    Xti = Yj(:,cIdx);

    % training network
    [network, trainInfo] = trainNetwork(Xti.', Xt.', net.layers, options);
    time = toc(ticH);

    trainedNet.nodeNum = nodeNum;
    trainedNet.exNum = exNum;
    trainedNet.network = network;
    trainedNet.trainInfo = trainInfo;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole nVARNN network! t = ' num2str(time) 's']);
end