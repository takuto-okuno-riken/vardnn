%%
% Caluclate Convergent Cross Mapping
% returns CCM causality index (CCM).
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  E            embedding dimension (default:3)
%  tau          time delay used in the phase-space reconstruction (default:1)
%  L            library of L points from the shadow manifolds (default:[]);
%  sampling     'linear' or 'random' (default:'linear')
%  isFullNode   return both node & exogenous causality matrix (optional)

% Before using this function, download matlab CCM codes from
% https://jp.mathworks.com/matlabcentral/fileexchange/52964-convergent-cross-mapping
% and add a path "CCM_L_M" folder.

function [CCM] = calcConvCrossMap__(X, exSignal, nodeControl, exControl, E, tau, isFullNode)
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

    % calc CCM
    CCM = nan(nodeNum, nodeMax);
    for i=1:nodeNum
        for j=i:nodeMax
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            [ SugiCorr , SugiR , LM , SugiY , SugiX , origY , origX ] = SugiLM(Y(i,:).', Y(j,:).', tau, E);
            CCM(i,j) = SugiCorr(1);
            if i~=j && j<=nodeNum && any(nodeControl(i,j,:),'all')
                CCM(j,i) = SugiCorr(2);
            end
        end
    end
end
