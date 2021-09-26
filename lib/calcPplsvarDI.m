%%
% Caluclate pPLSVAR (pairwise PLS Vector Auto-Regression) DI
% returns pPLSVAR DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          pPLSVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub, coeff] = calcPplsvarDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    inputNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + net.exNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc pairwise PCVAR
    DI = nan(nodeNum, nodeMax);
    coeff = nan(nodeNum, nodeMax);
    DIsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            b = net.bvec{i,j};
            z = sum(b);
            DIsub(i,j,1) = z;
            [~,idx] = find(control(i,j,:)==1);

            zj = z;
            for k=1:lags
                s = j + (k-1)*inputNum;
                bIdx = find(idx==s);
                if ~isempty(bIdx)
                    zj = zj - b(1+bIdx);
                end
            end
            DIsub(i,j,2) = z - zj;
            coeff(i,j) = zj; % actually this is sum of b(p+2:end)
            DI(i,j) = abs(coeff(i,j));
        end
    end
end

