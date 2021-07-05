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
    nodeInNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end

    % calc mPCVAR DI
    DI = nan(nodeNum,nodeMax);
    coeff = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeInNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:lags
            idx = [idx, nodeIdx+nodeInNum*(k-1), exIdx+nodeInNum*(k-1)];
        end
        idxList = [nodeIdx, exIdx];
        nlen = length(idxList);

        % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        mc = net.maxComp{i};
        mu = net.mu{i};

        Xti = ones(1,length(idx));
        score = (Xti - mu) / net.coeff{i}.';
        subScore = [score(:,1:mc), 1];
        z = subScore * net.bvec{i};
        DIsub(i,1) = z;

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            Xtj = Xti;
            bIdx = find(idxList==j);
            for k=1:lags
                Xtj(bIdx+nlen*(k-1)) = 0;
            end
            scorej = (Xtj - mu) / net.coeff{i}.';
            subScorej = [scorej(:,1:mc), 1];
            zj = subScorej * net.bvec{i};

            DI(i,j) = abs(z - zj);
            coeff(i,j) = z - zj;
            DIsub(i,j+1) = zj;
        end
    end
end

