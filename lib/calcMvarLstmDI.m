%%
% Calculate multivariate VAR LSTM Directional Influence matrix (DI) and impaired node signals (DIsub)
% returns mVAR LSTM DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained multivariate VAR LSTM network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = calcMvarLstmDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end

    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    if isfield(net, 'exNum'), exNum = net.exNum; else exNum = size(net.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2) - nodeNum; end % for compatibility
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc multivariate VAR LSTM DI
    DI = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    for i=1:nodeNum
        Xi = ones(nodeNum + exNum, lags);
        iNodeControl = squeeze(nodeControl(i,:,:));
        iExControl = squeeze(exControl(i,:,:));
        if size(nodeControl,3)<=1, iNodeControl = iNodeControl.'; end
        if size(exControl,3)<=1, iExControl = iExControl.'; end
        Xi(1:nodeNum,:) = Xi(1:nodeNum,:) .* iNodeControl;
        Xi(nodeNum+1:end,:) = Xi(nodeNum+1:end,:) .* iExControl;

        % predict 
        DIsub(i,1)  = predict(net.nodeNetwork{i}, Xi, 'ExecutionEnvironment', 'cpu');

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            Xtj = Xi;
            Xtj(j,:) = 0;

            % predict 
            DIsub(i,j+1) = predict(net.nodeNetwork{i}, Xtj, 'ExecutionEnvironment', 'cpu');
            DI(i,j) = abs(DIsub(i,1) - DIsub(i,j+1));
        end
    end
end

