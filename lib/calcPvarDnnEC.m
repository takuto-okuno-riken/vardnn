%%
% Caluclate pairwise VAR DNN-EC
% returns pairwise VAR DNN EC matrix (EC) and impaired node signals (ECsub)
% input:
%  net          trained Pairwised VAR DNN network structure
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [EC, ECsub] = calcPvarDnnEC(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    nodeInNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % calc Pairwised DNN granger causality
    EC = nan(nodeNum, nodeMax);
    ECsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            nodeInput = ones(2*lags,1);

            % predict 
            ECsub(i,j,1) = predict(net.nodeNetwork{i,j}, nodeInput);
            
            % imparement node signals
            impInput = nodeInput;
            impInput(lags+1:end,:) = 0;

            % predict 
            ECsub(i,j,2) = predict(net.nodeNetwork{i,j}, impInput);
            EC(i,j) = abs(ECsub(i,j,1)-ECsub(i,j,2));
        end
    end
end

