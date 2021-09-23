%%
% Caluclate mPCVAR (multivaliate Principal Component Vector Auto-Regression) Mean Inpact Value (MIV)
% returns Mean Inpact Value matrix (MIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPCVAR network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function MIV = calcMpcvarMIV(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    nodeInNum = nodeNum + net.exNum;
    p = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end

    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-p, p*nodeInNum);
    for k=1:p
        Yj(:,1+nodeInNum*(k-1):nodeInNum*k) = Y(1+k:sigLen-p+k,:);
    end

    MIV = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        nodeIdx = [1:nodeNum];
        if ~isempty(nodeControl)
            [~,nodeIdx] = find(nodeControl(i,:)==1);
        end
        exIdx = [nodeNum+1:nodeInNum];
        if ~isempty(exControl)
            [~,exIdx] = find(exControl(i,:)==1);
            exIdx = exIdx + nodeNum;
        end
        idx = [];
        for k=1:p
            idx = [idx, nodeIdx+nodeInNum*(k-1), exIdx+nodeInNum*(k-1)];
        end
        idxList = [nodeIdx, exIdx];
        nlen = length(idxList);
        Xti = Yj(:,idx);

        % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        mc = net.maxComp{i};
        mu = net.mu{i};

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            Xtj1 = Xti;
            Xtj2 = Xti;
            bIdx = find(idxList==j);
            for k=1:p
                Xtj1(:,bIdx+nlen*(k-1)) = Xtj1(:,bIdx+nlen*(k-1)) * 1.1;
                Xtj2(:,bIdx+nlen*(k-1)) = Xtj2(:,bIdx+nlen*(k-1)) * 0.9;
            end
            scorej1 = (Xtj1 - mu) / net.coeff{i}.';
            scorej2 = (Xtj2 - mu) / net.coeff{i}.';
            pcXti1 = [scorej1(:,1:mc), ones(sigLen-p,1)]; % might not be good to add bias
            pcXti2 = [scorej2(:,1:mc), ones(sigLen-p,1)]; % might not be good to add bias
            IV1 = pcXti1 * net.bvec{i};
            IV2 = pcXti2 * net.bvec{i};

            MIV(i,j) = mean(IV1 - IV2);
        end
    end
end

