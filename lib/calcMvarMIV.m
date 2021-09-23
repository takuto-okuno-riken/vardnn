%%
% Caluclate mVAR (multivaliate Vector Auto-Regression) Mean Inpact Value (MIV)
% returns mVAR MIV matrix (MIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained multivariate VAR network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function MIV = calcMvarMIV(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);

    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;

    p = net.lags;
    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-p, p*nodeMax);
    for k=1:p
        Yj(:,1+nodeMax*(k-1):nodeMax*k) = Y(1+k:sigLen-p+k,:);
    end

    % calc mVAR DI
    MIV = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeNum+exNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:p
            idx = [idx, nodeIdx+nodeMax*(k-1), exIdx+nodeMax*(k-1)];
        end

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            Yj1 = Yj; Yj2 = Yj;
            for k=1:p
                Yj1(:,j+nodeMax*(k-1)) = Yj(:,j+nodeMax*(k-1)) * 1.1;
                Yj2(:,j+nodeMax*(k-1)) = Yj(:,j+nodeMax*(k-1)) * 0.9;
            end
            Xti1 = [Yj1(:,idx), ones(sigLen-p,1)];
            Xti2 = [Yj2(:,idx), ones(sigLen-p,1)];

            IV1 = Xti1 * net.bvec{i};
            IV2 = Xti2 * net.bvec{i};
            MIV(i,j) = mean(IV1 - IV2);
        end
    end
end

