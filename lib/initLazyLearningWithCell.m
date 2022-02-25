%%
% init Lazy Learning structure
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for kNNS (default:3)
%  uniqueDecimal   taking unique value from conjuncted time series (option)

function LL = initLazyLearningWithCell(CX, CexSignal, nodeControl, exControl, lags, uniqueDecimal)
    if nargin < 6, uniqueDecimal = 0; end
    if nargin < 5, lags = 3; end
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

    % set node input
    for i=1:cxNum
        % set node input
        CY{i} = single(CX{i});
        if exNum > 0
            CY{i} = [CY{i}; single(CexSignal{i})];
        end
        CY{i} = flipud(CY{i}.'); % need to flip signal
    end

    LL.nodeNum = nodeNum;
    LL.exNum = exNum;
    LL.lags = lags;
    LL.Y = [];
    LL.CY = CY;
    LL.uniqueDecimal = uniqueDecimal;
end
