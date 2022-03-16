%%
% Return trained multivariate VAR DNN network
% input :
%  CX            cells of multivariate time series matrix {node x time series}
%  CexSignal     cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input) (optional)
%  net           mVAR DNN network structure
%  options       training options

function trainedNet = trainMvarDnnNetworkWithCell(CX, CexSignal, nodeControl, exControl, net, options)
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    sigLen = size(CX{1},2);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    lags = net.lags;
    inputNum = nodeNum + exNum;
    uniqueDecimal = net.uniqueDecimal;

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % training whole multivariate VAR DNN network
    disp('start training whole multivariate VAR DNN network');
    ticH = tic;

    % set vector auto-regression (VAR) inputs
    idxs = cell(nodeNum,1);
    for n=1:nodeNum
        [~,idxs{n}] = find(control(n,:,:)==1);
    end
    allInLen = 0;
    for i=1:cxNum
        allInLen = allInLen + size(CX{i},2) - lags;
    end

    % calculate mean and covariance of each node
    Y = [];
    for i=1:cxNum, Y = [Y, CX{i}]; end
    cxM = mean(Y.');
    cxCov = cov(Y.');
    Y = []; % memory clear

    nodeLayers = net.nodeLayers;
    nodeNetwork = cell(nodeNum,1);
    trainInfo = cell(nodeNum,1);
    rvec = cell(nodeNum,1);

%    for n=1:nodeNum
    parfor n=1:nodeNum    % for parallel processing
        if isempty(nodeLayers{n}), continue; end
        disp(['training node ' num2str(n)]);
        Xt = single(nan(allInLen,1));
        Xti = single(nan(allInLen,length(idxs{n})));
        xts = 1;

        % this is redundant for every node. but it is necessary to avoid
        % too much memory consumption
        for i=1:cxNum
            % set node input
            Y = single(CX{i});
            if exNum > 0
                Y = [Y; single(CexSignal{i})];
            end
            Y = flipud(Y.'); % need to flip signal

            sLen = size(Y,1);
            sl = sLen-lags;
            Yj = single(zeros(sl, lags*inputNum));
            for p=1:lags
                Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sl+p,:);
            end
            Xt(xts:xts+sl-1,:) = Y(1:sl,n);
            Xti(xts:xts+sl-1,:) = Yj(:,idxs{n});
            xts = xts + sl;
        end
        Y = [];  % clear memory
        Yj = []; % clear memory

        if uniqueDecimal > 0
            A = int32([Xt, Xti] / uniqueDecimal);
            A = unique(A,'rows');
            perm = randperm(size(A,1)); % need to shuffle for DNN training
            Xt = single(A(perm,1)) * uniqueDecimal;
            Xti = single(A(perm,2:end)) * uniqueDecimal;
            A = []; % clear memory
        end

        % training network
        [nodeNetwork{n}, trainInfo{n}] = trainNetwork(Xti.', Xt.', nodeLayers{n}, options);
        
        % calculate residuals for surrogate
        A = predict(nodeNetwork{n}, Xti.', 'ExecutionEnvironment', 'cpu');
        rvec{n} = Xt.' - A;
        A = []; % clear memory
    end
    time = toc(ticH);
    trainedNet = net;
    trainedNet.nodeNum = nodeNum; % for compatibility
    trainedNet.exNum = exNum; % for compatibility
    trainedNet.sigLen = sigLen;
    trainedNet.cxM = cxM;
    trainedNet.cxCov = cxCov;
    trainedNet.nodeNetwork = nodeNetwork;
    trainedNet.trainInfo = trainInfo;
    trainedNet.trainTime = time;
    trainedNet.trainOptions = options;
    trainedNet.rvec = rvec;
    disp(['finish training whole multivariate VAR DNN network! t = ' num2str(time) 's']);
end