%%
% Estimate hidden neurons and initial weight and create multivariate VAR DNN
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:1)
%  uniqueDecimal   taking unique value from conjuncted time series (option)
%  activateFunc    activation function for each layer (default:@reluLayer)
%  initWeightFunc  initializing weight function (default:[])
%  initWeightParam parameters for initializing weight function (default:[])
%  initBias        initializing bias value (default:0)
%             For uniform distribution, bias = 0 and empty initial weight is better
%             For fMRI BOLD signal, bias = 0.5 and rough initial weight is better

function net = initMvarDnnNetworkWithCell(CX, CexSignal, nodeControl, exControl, lags, uniqueDecimal, activateFunc, initWeightFunc, initWeightParam, initBias)
    if nargin < 10, initBias = 0; end
    if nargin < 9, initWeightParam = []; end
    if nargin < 8, initWeightFunc = []; end
    if nargin < 7, activateFunc = @reluLayer; end
    if nargin < 6, uniqueDecimal = 0; end
    if nargin < 5, lags = 1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, CexSignal = {}; end
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    inputNum = nodeNum + exNum;

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    estimateInNum = ceil((sum(nodeControl,'all') + sum(exControl,'all')) / nodeNum);

    % set vector auto-regression (VAR) inputs
    idxs = cell(nodeNum,1);
    for n=1:nodeNum
        [~,idxs{n}] = find(control(n,:,:)==1);
    end
    allInLen = 0;
    for i=1:cxNum
        allInLen = allInLen + size(CX{i},2) - lags;
    end

    sigLens = nan(nodeNum,1);
    for n=1:4:nodeNum
%    parfor n=1:4:nodeNum
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
            Xt = single(A(:,1)) * uniqueDecimal;
            A = []; % clear memory
        end

        sigLens(n) = size(Xt,1);
        Xt = [];  % clear memory
        Xti = []; % clear memory
    end
    
    % estimate neuron number of hidden layers
    hiddenNums = estimateHiddenNeurons(estimateInNum, nanmean(sigLens));
    
    % layer parameters
    net = createMvarDnnNetwork(nodeNum, exNum, hiddenNums, lags, nodeControl, exControl, activateFunc, initWeightFunc, initWeightParam, initBias);
    net.uniqueDecimal = uniqueDecimal;
end
