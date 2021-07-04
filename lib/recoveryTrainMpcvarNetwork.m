%%
% Return trained multivariate PCVAR network
% input :
%  X             multivariate time series matrix (node x time series)
%  exSignal      multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl   node control matrix (node x node) (optional)
%  exControl     exogenous input control matrix for each node (node x exogenous input)
%  net           network structure
%  maeMin        minimum MAE of re-training (default:0.008)
%  retrainRange  error range rate(0-1) of re-training (default:[0.1, 0.3])
%  retrainMax    muximum re-training number (default:10)

function [trainedNet, time, mae, errs] = recoveryTrainMpcvarNetwork(X, exSignal, nodeControl, exControl, net, maeMin, retrainRange, retrainMax)
    if nargin < 8, retrainMax = 10; end
    if nargin < 7, retrainRange = [0.1, 0.3]; end
    if nargin < 6, maeMin = 0.008; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    expTh = net.explainedTh * 100;
    trainedNet = net;

    ticH = tic; % start stop watch

    % simulate network with 1st frame & exogenous input signal
    [S, time] = simulateMpcvarNetwork(X, exSignal, nodeControl, exControl, net);
    [mae, ~, errs] = getTwoSignalsError(X, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);

    errX = mean(abs(errs),1);
    errMax = max(errX);
    errRange(1) = retrainRange(1) * errMax;
    errRange(2) = retrainRange(2) * errMax;
%    fg = figure; plot(errX);

    errIdx = find(errRange(1) < errX & errX < errRange(2));
    if mae < maeMin
        time = toc(ticH); return;
    end
    if isempty(errIdx) % bad case
        [B, I] = sort(errX);
        for i=1:length(B)
            if B(i)>0, errIdx = [errIdx, I(i)]; end
            if length(errIdx)>2, break; end
        end
    end
    lastmae = mae;

    % set input control indexes
    nodeMax = nodeNum + exNum;
    idxs = {};
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeNum+exNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:lags
            idx = [idx, nodeIdx+nodeMax*(k-1), exIdx+nodeMax*(k-1)];
        end
        idxs{i} = idx;
    end

    % recovery training whole DLCM network
    Y = [X; exSignal];
    baseS = [];
    inputHist = [];
    teachHist = [];
    for k=1:lags, baseS = [baseS; Y(:,lags-(k-1):end-k)]; end
    baseTeach = X(:,lags+1:end);
    lastY = [S; exSignal];

    disp('start recovery training whole network');
    for k=1:retrainMax
        if errIdx(end)==sigLen, errIdx(end)=[]; end
        errInput = [];
        for j=1:length(errIdx)
            eY = [];
            for i=1:lags, eY = [eY; lastY(:,errIdx(j)-(i-1));]; end
            errInput = [errInput, eY];
        end
        errTeach = X(:,errIdx+1);
        nodeInput = [baseS, inputHist, errInput];
        nodeTeach = [baseTeach, teachHist, errTeach];

        ticH2 = tic;
        coeff = cell(nodeNum,1);
        score = cell(nodeNum,1);
        latent = cell(nodeNum,1);
        explained = cell(nodeNum,1);
        mu = cell(nodeNum,1);
        maxComp = cell(nodeNum,1);
        b = cell(nodeNum,1);
        bint = cell(nodeNum,1);
        r = cell(nodeNum,1);
        rint = cell(nodeNum,1);
        stats = cell(nodeNum,1);

        for i=1:nodeNum
            % vector auto-regression (VAR)
            Xt = nodeTeach(i,:).';
            Xti = nodeInput(idxs{i},:).';

            % apply the Principal Component Regress function
            [coeff{i},score{i},latent{i},~,explained{i},mu{i}] = pca(Xti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

            % find 99% component range
            expTotal = 0;
            maxComp = size(score,2);
            for j=1:size(Xti,2)
                expTotal = expTotal + explained{i}(j);
                if expTotal >= expTh
                    maxComp{i} = j;
                    break;
                end
            end
            pcXti = [score{i}(:,1:maxComp{i}), ones(size(Xt,1),1)]; % might not be good to add bias
            [b{i},bint{i},r{i},rint{i},stats{i}] = regress(Xt, pcXti);
        end
        net.coeff = coeff;
        net.score = score;
        net.latent = latent;
        net.explained = explained;
        net.mu = mu;
        net.maxComp = maxComp;
        net.bvec = b;
        net.bint = bint;
        net.rvec = r;
        net.rint = rint;
        net.stats = stats;
        time = toc(ticH2);
        disp(['recovery training (' num2str(k) '): train time=' num2str(time)]);

        % simulate again network with 1st frame & exogenous input signal
        [S, time] = simulateMpcvarNetwork(X, exSignal, nodeControl, exControl, net);
        [mae, ~, errs2] = getTwoSignalsError(X, S);
        disp(['recovery training (' num2str(k) '): simulation time=' num2str(time) ', mae=' num2str(mae)]);

        errX = mean(abs(errs2),1);
        errs = cat(3, errs, errs2);
        errIdx = find(errRange(1) < errX & errX < errRange(2));
        if lastmae > mae
            lastmae = mae;
            lastY = [S; exSignal];
            inputHist = [inputHist, errInput];
            teachHist = [teachHist, errTeach];
            trainedNet.coeff = coeff;
            trainedNet.score = score;
            trainedNet.latent = latent;
            trainedNet.explained = explained;
            trainedNet.mu = mu;
            trainedNet.maxComp = maxComp;
            trainedNet.bvec = b;
            trainedNet.bint = bint;
            trainedNet.rvec = r;
            trainedNet.rint = rint;
            trainedNet.stats = stats;
        end
        if mae < maeMin
            break;
        end
        if isempty(errIdx) % bad case
            [B, I] = sort(errX);
            for i=1:length(B)
                if B(i)>0, errIdx = [errIdx, I(i)]; end
                if length(errIdx)>2, break; end
            end
        end
    end
    time = toc(ticH);
    trainedNet.recoverTrainTime = time;
    disp(['finish recovery training whole network! t = ' num2str(time) 's']);
end