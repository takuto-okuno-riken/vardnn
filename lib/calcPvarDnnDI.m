%%
% Calculate pairwise VAR DNN-DI
% returns pairwise VAR DNN DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained Pairwise VAR DNN network structure
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [DI, DIsub] = calcPvarDnnDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    nodeNum = net.nodeNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + net.exNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc Pairwise DNN DI
    DI = nan(nodeNum, nodeMax);
    DIsub = nan(nodeNum, nodeMax, 2);
    for i=1:nodeNum
        [~,idx] = find(control(i,i,:)==1);
        Yi = ones(1,lags);
        Xti = Yi(:,idx);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            [~,idx] = find(control(i,j,:)==1);
            Yj = ones(1,lags);
            Xtj = Yj(:,idx);

            % predict 
            DIsub(i,j,1) = predict(net.nodeNetwork{i,j}, [Xti,Xtj].');
            
            % imparement node signals
            Yj(:) = 0;
            Xtj = Yj(:,idx);

            % predict 
            DIsub(i,j,2) = predict(net.nodeNetwork{i,j}, [Xti,Xtj].');
            DI(i,j) = abs(DIsub(i,j,1)-DIsub(i,j,2));
        end
    end
end

