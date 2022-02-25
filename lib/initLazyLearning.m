%%
% init Lazy Learning structure
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for kNNS (default:3)

function LL = initLazyLearning(X, exSignal, nodeControl, exControl, lags)
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    exNum = size(exSignal,1);

    % set node input
    Y = [X; exSignal];    
    Y = flipud(Y.'); % need to flip signal

    LL.nodeNum = nodeNum;
    LL.exNum = exNum;
    LL.lags = lags;
    LL.Y = Y;
end