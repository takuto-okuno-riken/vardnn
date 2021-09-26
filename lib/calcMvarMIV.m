%%
% Caluclate mVAR (multivaliate Vector Auto-Regression) Mean Impact Value (MIV)
% returns mVAR MIV matrix (MIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained multivariate VAR network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [MIV, MAIV] = calcMvarMIV(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeNum + exNum; end

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end

    % calc mVAR MIV
    MIV = nan(nodeNum,nodeMax);
    MAIV = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end
            Yj1 = Yj; Yj2 = Yj;
            for k=1:lags
                Yj1(:,j+nodeMax*(k-1)) = Yj(:,j+nodeMax*(k-1)) * 1.1;
                Yj2(:,j+nodeMax*(k-1)) = Yj(:,j+nodeMax*(k-1)) * 0.9;
            end
            Xti1 = [Yj1(:,idx), ones(sigLen-lags,1)];
            Xti2 = [Yj2(:,idx), ones(sigLen-lags,1)];

            IV1 = Xti1 * net.bvec{i};
            IV2 = Xti2 * net.bvec{i};
            MIV(i,j) = mean(IV1 - IV2);
            MAIV(i,j) = mean(abs(IV1 - IV2));
        end
    end
end

