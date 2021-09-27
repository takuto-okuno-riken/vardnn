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
    inputNum = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal
    
    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    nodeLayers = cell(nodeNum,1);
    coeff = cell(nodeNum,1);
    score = cell(nodeNum,1);
    latent = cell(nodeNum,1);
    explained = cell(nodeNum,1);
    mu = cell(nodeNum,1);

    % set input signals
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end
    for i=1:nodeNum
%    parfor i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);

        % PCA inputs and VARDNN teaching signals
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);
        
        % apply the Principal Component Regress function
        [coeff{i},score{i},latent{i},~,explained{i},mu{i}] = pca(Xti); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

        % estimate neuron number of hidden layers
        pcNodeNum = size(score{i},2);
        hiddenNums = estimateHiddenNeurons(pcNodeNum, sigLen);

        % init DNN layers
        nodeLayers{i} = createMvarDnnLayers(hiddenNums, ones(1,pcNodeNum), [], activateFunc, initWeightFunc, initWeightParam, initBias, i);
    end
    net.version = 1.1;
    net.nodeNum = nodeNum;
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