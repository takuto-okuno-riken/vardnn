%%
% Caluclate Convergent Cross Mapping - FC (subtract FC)
% returns CCM causality index (CCM), FC p-values (Pfc) and CCM p-values (Pccm).
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  E            embedding dimension (default:3)
%  tau          time delay used in the phase-space reconstruction (default:1)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [CCM, Pfc, Pccm] = calcConvCrossMapSubFC(X, exSignal, nodeControl, exControl, E, tau, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, tau = 1; end
    if nargin < 5, E = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, E);

    [CCM, Pccm] = calcConvCrossMap_(X, exSignal, nodeControl, exControl, E, tau, isFullNode);
    [FC, Pfc] = corr(Y.', Y.');

    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end
            
            CCM(i,j) = CCM(i,j) - FC(i,j);
        end
    end
end
