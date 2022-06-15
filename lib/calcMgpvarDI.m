%%
% Calculate mGPVAR (multivaliate Gaussian Processes Vector Auto-Regression) DI
% returns mGPVAR DI matrix (DI), impaired node signals (DIsub) and regression
% coefficient matrix (coeff)
% input:
%  net          mGPVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [DI, DIsub] = calcMgpvarDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    exNum = net.exNum;
    inputNum = nodeNum + exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc multivariate SVM VAR DI
    DI = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    nodeInputOrg = ones((nodeNum + exNum)*lags, 1);
    for i=1:nodeNum
        if isempty(net.mdl{i}), continue; end
        [~,idx] = find(control(i,:,:)==1);
        nodeInput = nodeInputOrg(idx,1);

        % predict 
        DIsub(i,1) = predict(net.mdl{i}, nodeInput.');

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            nodeInput = nodeInputOrg;
            for p=1:lags, nodeInput(j+(nodeNum + exNum)*(p-1),:) = 0; end
            nodeInput = nodeInput(idx,1);

            % predict 
            DIsub(i,j+1) = predict(net.mdl{i}, nodeInput.');
            DI(i,j) = abs(DIsub(i,1) - DIsub(i,j+1));
        end
    end
end

