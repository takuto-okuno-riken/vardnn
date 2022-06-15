%%
% Calculate multivariate VAR DNN Directional Influence 3D matrix (DI) and impaired node signals (DIsub)
% returns mVAR DNN DI 3D matrix (DI) and impaired node signals (DIsub)
% input:
%  net          trained multivariate VAR DNN network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [DI, DIsub] = calcMvarDnnDI3(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end

    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    if isfield(net, 'exNum'), exNum = net.exNum; else exNum = size(net.nodeNetwork{1, 1}.Layers(2, 1).Weights, 2) - nodeNum; end % for compatibility
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isfield(net, 'version'), version = net.version; else version = 1; end
    inputNum = nodeNum + exNum;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc multivariate VAR DNN DI
    DI = nan(nodeNum,nodeMax,lags);
    DIsub = nan(nodeNum,nodeMax+1,lags);
    nodeInput = ones(inputNum*lags, 1);
    for i=1:nodeNum
        if isempty(net.nodeNetwork{i}), continue; end
        [~,idx] = find(control(i,:,:)==1);
        Xti = nodeInput(idx,1);

        % predict 
        DIsub(i,1,1) = predict(net.nodeNetwork{i}, Xti, 'ExecutionEnvironment', 'cpu');

        % imparement node signals
        for j=1:nodeMax
            for p=1:lags
                Xtp = nodeInput;
                Xtp(j+inputNum*(p-1),:) = 0;
                Xtp = Xtp(idx,1);

                % predict 
                DIsub(i,j+1,p) = predict(net.nodeNetwork{i}, Xtp, 'ExecutionEnvironment', 'cpu');
                DI(i,j,p) = abs(DIsub(i,1,1) - DIsub(i,j+1,p));
            end
        end
    end
end

