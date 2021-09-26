%%
% Caluclate mPCVAR (multivaliate Principal Component Vector Auto-Regression) DI
% returns mPCVAR DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          mPCVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub, coeff] = calcMpcvarDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    inputNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc mPCVAR DI
    Yj = ones(1, lags*inputNum);
    DI = nan(nodeNum,nodeMax);
    coeff = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);
        Xti = Yj(:,idx);

        % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        mc = net.maxComp{i};
        mu = net.mu{i};

        score = (Xti - mu) / net.coeff{i}.';
        subScore = [score(:,1:mc), 1];
        z = subScore * net.bvec{i};
        DIsub(i,1) = z;

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            Ytj = Yj;
            for k=1:lags, Ytj(:,j+inputNum*(k-1)) = 0; end
            Xtj = Ytj(:,idx);

            scorej = (Xtj - mu) / net.coeff{i}.';
            subScorej = [scorej(:,1:mc), 1];
            zj = subScorej * net.bvec{i};

            DI(i,j) = abs(z - zj);
            coeff(i,j) = z - zj;
            DIsub(i,j+1) = zj;
        end
    end
end

