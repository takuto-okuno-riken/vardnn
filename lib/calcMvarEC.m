%%
% Caluclate mVAR (multivaliate Vector Auto-Regression) EC
% returns mVAR EC matrix (EC), impaired node signals (ECsub) and regression
% coefficient matrix (coeff)
% input:
%  net          mVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [EC, ECsub, coeff] = calcMvarEC(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % calc mVAR EC
    EC = nan(nodeNum,nodeMax);
    ECsub = nan(nodeNum,nodeMax+1);
    for i=1:nodeNum
        b = net.bvec{i};
        z = sum(b);
        ECsub(i,1) = z;

        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeNum+net.exNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idxList = [nodeIdx, exIdx];
        nlen = length(idxList);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            zj = z;
            bIdx = find(idxList==j);
            for k=1:lags
                zj = zj - b(bIdx+nlen*(k-1));
            end

            EC(i,j) = abs(z - zj); % actually this is sum of b(bIdx+nlen*(k-1))
            ECsub(i,j+1) = zj;
        end
    end
    coeff = repmat(ECsub(:,1), [1 size(EC,2)]) - ECsub(:,2:end);
end

