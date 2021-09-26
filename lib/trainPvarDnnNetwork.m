%%
% Return trained pairwise VAR DNN network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           Pairwise VAR DNN network structure
%  options       training options

function trainedNet = trainPvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options)
    global dnnInitWeights;
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    sigLen = size(X,2);
    trainedNet = net;
    lags = net.lags;

    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    Y = flipud(Y.'); % need to flip signal

    % training whole Pairwise VAR DNN network
    disp('start training whole pairwise VAR DNN network');
    ticH = tic;
    nodeLayers = trainedNet.nodeLayers;
    nodeNetwork = cell(nodeNum,nodeMax);
    trainInfo = cell(nodeNum,nodeMax);
    initWeights = cell(nodeNum,nodeMax);
    for i=1:nodeNum
        [~,idx] = find(control(i,i,:)==1);
        Xt = Y(1:sigLen-lags,i);
        Yi = zeros(sigLen-lags, lags);
        for k=1:lags, Yi(:,k) = Y(1+k:sigLen-lags+k,i); end
        Xti = Yi(:,idx);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum &&  nodeControl(i,j,1) == 0, continue; end
            if j>nodeNum &&  exControl(i,j-nodeNum,1) == 0, continue; end

            disp(['training node ' num2str(i) '-' num2str(j)]);
            [~,idx] = find(control(i,j,:)==1);
            Yj = zeros(sigLen-lags, lags);
            for k=1:lags, Yj(:,k) = Y(1+k:sigLen-lags+k,j); end
            Xtj = Yj(:,idx);
            [nodeNetwork{i,j}, trainInfo{i,j}] = trainNetwork([Xti,Xtj].', Xt.', nodeLayers{i,j}, options);
            initWeights{i,j} = dnnInitWeights;
        end
    end
    time = toc(ticH);
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.initWeights = initWeights;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    disp(['finish training whole Pairwise VAR DNN network! t = ' num2str(time) 's']);
end