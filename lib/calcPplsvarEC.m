%%
% Caluclate pPLSVAR (pairwise PLS Vector Auto-Regression) EC
% returns pPLSVAR EC matrix (EC) and impaired node signals (ECsub)
% input:
%  net          pPLSVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [EC, ECsub] = calcPplsvarEC(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;
    p = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % calc pairwise PCVAR
    EC = nan(nodeNum, nodeMax);
    ECsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            ECsub(i,j,1) = sum(net.bvec{i,j});
            ECsub(i,j,2) = sum(net.bvec{i,j}(1:p+1));
            EC(i,j) = abs(ECsub(i,j,1)-ECsub(i,j,2)); % actually this is sum of b(p+2:end)
        end
    end
end

