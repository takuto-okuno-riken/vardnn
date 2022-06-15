%%
% Calculate Nonlinear VAR Neural Network (NVARNN) Mean Impact Value of 3D matrix (Full Input)
% returns Nonlinear VAR Neural Network (NVARNN) Mean Impact Value 3D matrix (MIV) and Mean Absolute Impact Value 3D matrix (MAIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (1 x node) (optional)
%  exControl    exogenous input control matrix for each node (1 x exogenous input) (optional)
%  net          trained nVARNN network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [MIV, MAIV] = calcNvarnnMIV3(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end

    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    sigLen = size(X,2);
    inputNum = nodeNum + exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (1 x node x lags)
    if isempty(nodeControl)
        nodeControl = ones(1,nodeNum,lags);
    end
    if isempty(exControl)
        exControl = ones(1,exNum,lags);
    end
    [~, ~, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    cIdx = find(control(:)==1);

    % set training data
    Yj = zeros(sigLen-lags, lags*inputNum);
    for p=1:lags
        Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sigLen-lags+p,:);
    end

    % calc nVARNN MIV
    MIV = nan(nodeNum, nodeMax, lags);
    MAIV = nan(nodeNum, nodeMax, lags);

    % imparement node signals
    for j=1:nodeMax
        for p=1:lags
            Xj1 = Yj;
            Xj2 = Yj;
            n = j + inputNum * (p-1);
            Xj1(:,n) = Xj1(:,n) * 1.1; 
            Xj2(:,n) = Xj2(:,n) * 0.9; 
            Xj1 = Xj1(:,cIdx);
            Xj2 = Xj2(:,cIdx);

            % predict 
            N1 = predict(net.network, Xj1.', 'ExecutionEnvironment', 'cpu');
            N2 = predict(net.network, Xj2.', 'ExecutionEnvironment', 'cpu');
            MIV(:,j,p) = mean(N1 - N2,2);
            MAIV(:,j,p) = mean(abs(N1 - N2),2);
        end
    end
end

