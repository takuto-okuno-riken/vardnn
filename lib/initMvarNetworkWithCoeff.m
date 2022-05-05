%%
% create mVARX network with coefficient
% input:
%  coefficient     VAR coefficient (node x [(node, ex) x lags, 1])
%  nodeNum         node number
%  exNum           exogenous node number
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)

function net = initMvarNetworkWithCoeff(coefficient, nodeNum, exNum, nodeControl, exControl, lags)
    if nargin < 6, lags = 3; end
    if nargin < 5, exControl = []; end
    if nargin < 4, nodeControl = []; end

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    b = cell(nodeNum,1);
    r = cell(nodeNum,1);
    T = cell(nodeNum,1);

    % first, calculate vector auto-regression (VAR) without target
    for i=1:nodeNum
        b{i} = coefficient(i,:).';
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.bvec = b;
    net.rvec = r;
    net.Tvec = T;
end