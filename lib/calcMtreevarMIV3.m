%%
% Caluclate multivariate Tree VAR Mean Impact Value of 3D matrix (Full Input)
% returns multivariate Tree VAR Mean Impact Value 3D matrix (MIV) and Mean Absolute Impact Value 3D matrix (MAIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained multivariate Tree VAR network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [MIV, MAIV] = calcMtreevarMIV3(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end

    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    sigLen = size(X,2);
    inputNum = nodeNum + exNum;
    if isfield(net, 'lags'), lags = net.lags; else lags = 1; end
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    Y = flipud(Y.'); % need to flip signal

    % set training data
    Yj = zeros(sigLen-lags, lags*inputNum);
    for p=1:lags
        Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sigLen-lags+p,:);
    end

    % calc mTree VAR MIV
    MIV = nan(nodeNum, nodeMax, lags);
    MAIV = nan(nodeNum, nodeMax, lags);
    for i=1:nodeNum
        if isempty(net.trees{i}), continue; end
        [~,idx] = find(control(i,:,:)==1);
        
        % imparement node signals
        for j=1:nodeMax
            for p=1:lags
                Xj1 = Yj;
                Xj2 = Yj;
                n = j + inputNum * (p-1);
                Xj1(:,n) = Xj1(:,n) * 1.1; 
                Xj2(:,n) = Xj2(:,n) * 0.9; 
                Xj1 = Xj1(:,idx);
                Xj2 = Xj2(:,idx);

                % predict 
                N1 = predict(net.trees{i}, Xj1);
                N2 = predict(net.trees{i}, Xj2);
                MIV(i,j,p) = mean(N1 - N2);
                MAIV(i,j,p) = mean(abs(N1 - N2));
            end
        end
    end
end

