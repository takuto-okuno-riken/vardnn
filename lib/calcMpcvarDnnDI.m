%%
% Calculate multivariate PC VAR DNN Directional Influence matrix (DI) and impaired node signals (DIsub)
% returns mPC VAR DNN DI matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained multivariate PC VAR DNN network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = calcMpcvarDnnDI(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end

    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    if isfield(net, 'exNum'), exNum = net.exNum; else exNum = size(net.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2) - nodeNum; end % for compatibility
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end

    % set input control indexes
    nodeNum = net.nodeNum;
    inputNum = nodeNum + net.exNum;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    % calc multivariate PC VAR DNN DI
    mu = net.mu;
    coeff = net.coeff;
    nodeNetwork = net.nodeNetwork;

    S1 = ones(inputNum*lags, 1);
    DI = nan(nodeNum,nodeMax);
    DIsub = nan(nodeNum,nodeMax+1);
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);
        S2 = S1(idx,1);
        Z = (S2.' - mu{i}) / coeff{i}.';
        % predict
        DIsub(i,1) = predict(nodeNetwork{i}, Z.', 'ExecutionEnvironment', 'cpu');

        % imparement node signals
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            S3 = S1;
            for p=1:lags, S3(j+inputNum*(p-1),:) = 0; end
            S3 = S3(idx,1);
            Z = (S3.' - mu{i}) / coeff{i}.';

            % predict 
            DIsub(i,j+1) = predict(nodeNetwork{i}, Z.', 'ExecutionEnvironment', 'cpu');
            DI(i,j) = abs(DIsub(i,1) - DIsub(i,j+1));
        end
    end
end

