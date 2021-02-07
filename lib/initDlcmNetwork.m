%%
% Estimate hidden neurons and initial weight and Create DLCM neural network
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

function net = initDlcmNetwork(X, exSignal, nodeControl, exControl, lags, activateFunc, initWeightFunc, initWeightParam, initBias)
    if nargin < 9, initBias = 0; end
    if nargin < 8, initWeightParam = []; end
    if nargin < 7, initWeightFunc = []; end
    if nargin < 6, activateFunc = @reluLayer; end
    if nargin < 5, lags = 1; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    net = initMvarDnnNetwork(X, exSignal, nodeControl, exControl, lags, activateFunc, initWeightFunc, initWeightParam, initBias);
end