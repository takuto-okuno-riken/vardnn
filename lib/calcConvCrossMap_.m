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
%  isFullNode   return both node & exogenous causality matrix (optional)

function [CCM] = calcConvCrossMap_(X, exSignal, nodeControl, exControl, E, tau, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, tau = 1; end
    if nargin < 5, E = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, E);

    % calc CCM
    Midx = cell(nodeMax,1);
    Mdist = cell(nodeMax,1);
    embtLen = sigLen - (E-1)*tau;
    Z = zeros(embtLen,E);
    for i=1:nodeMax
        for j=1:E
            Z(:,E-(j-1)) = Y(i,j:embtLen+(j-1));
        end
        % find K nearest neighbors of each Y time point on shadow manifold j
        [Midx{i}, Mdist{i}] = knnsearch(Z,Z,'K',E+2,'distance','euclidean');
    end
    
    CCM = nan(nodeNum, nodeMax);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % find K nearest neighbors on Yi from shadow manifold j
            MjIdx = Midx{j};
            Mjd = Mdist{j};
            Yi = Y(i,:);
            Yinn = Yi(MjIdx(:,2:E+2)+(E-1)*tau);

            % predict Yi feature points
            W = exp(-Mjd(:,2:E+2) ./ (Mjd(:,2) + 1e-8));
            W = W ./ sum(W, 2);
            Ypred = W .* Yinn;
            Ypred = sum(Ypred,2);

            % compare original Yi vs. predicted Yi
            % if it correlated well, j CCM cause i.
            Y1 = Y(i, 1+(E-1)*tau:end);
            CCM(i,j) = corr(Ypred, Y1');
        end
    end
end
