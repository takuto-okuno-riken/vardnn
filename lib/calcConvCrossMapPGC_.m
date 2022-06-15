%%
% Calculate Convergent Cross Mapping Pairwise Granger Causality
% returns CCM causality index (CCM) and p-values (P).
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  E            embedding dimension (default:3)
%  tau          time delay used in the phase-space reconstruction (default:1)
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [CCM, P] = calcConvCrossMapPGC_(X, exSignal, nodeControl, exControl, E, tau, alpha, isFullNode)
    if nargin < 8, isFullNode = 0; end
    if nargin < 7, alpha = 0.05; end
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
    stepAhead = (E-1)*tau;
    Knn = E+2;
%    if stepAhead>1, stepAhead=1; end % this affects a lot
    Midx = cell(nodeMax,1);
    Mdist = cell(nodeMax,1);
    embtLen = sigLen - (E-1)*tau;
    Z = zeros(embtLen,E);
    for i=1:nodeMax
        for j=1:E
            Z(:,E-(j-1)) = Y(i,j:embtLen+(j-1));
        end
        % find K nearest neighbors of each Y time point on shadow manifold j
        [Midx{i}, Mdist{i}] = knnsearch(Z,Z,'K',Knn,'distance','euclidean');
    end
    
    CCM = nan(nodeNum, nodeMax);
    P = nan(nodeNum, nodeMax);
    for i=1:nodeNum
        % find K nearest neighbors on Yi from shadow manifold i
        % and get N step ahead points
        MjIdx = Midx{i};
        Mjd = Mdist{i};
        Yi = Y(i,:);
        Aidx = MjIdx(:,2:Knn) + stepAhead;
        Yinn = Yi(Aidx);

        % predict Yi feature points
        W = exp(-Mjd(:,2:Knn) ./ (Mjd(:,2) + 1e-50));
        W = W ./ sum(W, 2);
        Ypred = W .* Yinn;
        YpredA = sum(Ypred,2);

        % apply the regress function and calc var of residuals
        Y1 = Y(i, 1+stepAhead:embtLen+stepAhead).';
        Xti = [YpredA ones(embtLen,1)];
        [b,bint,Xr] = regress(Y1,Xti);
        Vxt = var(Xr,1);
            
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % find K nearest neighbors on Yi from shadow manifold i&j
            % and get N step ahead points
            Z = zeros(embtLen,E*2);
            for k=1:E
                Z(:,E-(k-1)) = Y(i,k:embtLen+(k-1));
                Z(:,2*E-(k-1)) = Y(j,k:embtLen+(k-1));
            end
            [MjIdx, Mjd] = knnsearch(Z,Z,'K',Knn,'distance','euclidean');
            Aidx = MjIdx(:,2:Knn) + stepAhead;
            Yinn = Yi(Aidx);

            % predict Yi feature points
            W = exp(-Mjd(:,2:Knn) ./ (Mjd(:,2) + 1e-50));
            W = W ./ sum(W, 2);
            Ypred = W .* Yinn;
            YpredB = sum(Ypred,2);

            % apply the regress function and calc var of residuals
            Yti = [YpredB, ones(embtLen,1)];
            [b,bint,Yr] = regress(Y1,Yti);
            Vyt = var(Yr,1);
            if Vyt == 0
                 Vyt = 1e-15; % TODO: dummy to avoid inf return
            end
    
            CCM(i,j) = log(Vxt / Vyt);

            % TODO: calc F-statistic
            % https://en.wikipedia.org/wiki/F-test
            % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
            RSS1 = Xr'*Xr;  % p1 = p + 1
            RSS2 = Yr'*Yr;  % p2 = k
            k = E*2 + 1;
            F = ((RSS1 - RSS2)/E) / (RSS2 / (sigLen - k));
            P(i,j) = 1 - fcdf(F,E,(sigLen-k));
            cvFd = finv(1-alpha,E,(sigLen-k));
            h = F > cvFd;
        end
    end
end
