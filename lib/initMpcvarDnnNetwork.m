%%
% Estimate hidden neurons and initial weight and create multivariate PC VAR DNN
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:1)
%  activateFunc    activation function for each layer (default:@reluLayer)
%  initWeightFunc  initializing weight function (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = initMpcvarDnnNetwork(X, exSignal, nodeControl, exControl, lags, activateFunc, initWeightFunc, initWeightParam, initBias)
    if nargin < 9, initBias = []; end
    if nargin < 8, initWeightParam = []; end
    if nargin < 7, initWeightFunc = []; end
    if nargin < 6, activateFunc = @reluLayer; end
    if nargin < 5, lags = 1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % estimate PCA component number (with maximum node+exogenous)
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;
    [coeff,score,latent,tsquared,explained,mu] = pca(Y.');
    pcNodeNum = size(score,2);

    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(pcNodeNum, sigLen);
    
    % layer parameters
    p = lags;
    Y = flipud(Y.'); % need to flip signal

    nodeLayers = cell(nodeNum,1);
    coeff = cell(nodeNum,1);
    score = cell(nodeNum,1);
    latent = cell(nodeNum,1);
    explained = cell(nodeNum,1);
    mu = cell(nodeNum,1);

    % set input signals
    Yj = zeros(sigLen-p, p*nodeMax);
    for k=1:p
        Yj(:,1+nodeMax*(k-1):nodeMax*k) = Y(1+k:sigLen-p+k,:);
    end
    for i=1:nodeNum
%    parfor i=1:nodeNum
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
        for k=1:p
            idx = [idx, nodeIdx+nodeMax*(k-1), exIdx+nodeMax*(k-1)];
        end

        % PCA inputs and VARDNN teaching signals
        Xt = Y(1:sigLen-p,i);
        Xti = Yj(:,idx);
        
        % apply the Principal Component Regress function
        [coeff{i},score{i},latent{i},~,explained{i},mu{i}] = pca(Xti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

        % init DNN layers
        nodeLayers{i} = createMvarDnnLayers(hiddenNums, ones(1,size(score{i},2)), [], activateFunc, initWeightFunc, initWeightParam, initBias, i);
    end
    net.nodeNum = nodeNum;
    net.pcNodeNum = pcNodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.coeff = coeff;
    net.mu = mu;
%    net.score = score; % big size. calc this later
    net.nodeLayers = nodeLayers;
    net.nodeNetwork = cell(nodeNum,1);
    net.trainInfo = cell(nodeNum,1);
    net.initWeights = cell(nodeNum,1);
    net.trainTime = 0;
    net.recoverTrainTime = 0;
end