%%
% Caluclate pairwise VAR DNN-DI
% returns pairwise VAR DNN DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained Pairwised VAR DNN network structure
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [DI, DIsub] = calcPvarDnnDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % calc Pairwised DNN DI
    DI = nan(nodeNum, nodeMax);
    DIsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            nodeInput = ones(2*lags,1);

            % predict 
            DIsub(i,j,1) = predict(net.nodeNetwork{i,j}, nodeInput);
            
            % imparement node signals
            impInput = nodeInput;
            impInput(lags+1:end,:) = 0;

            % predict 
            DIsub(i,j,2) = predict(net.nodeNetwork{i,j}, impInput);
            DI(i,j) = abs(DIsub(i,j,1)-DIsub(i,j,2));
        end
    end
end

