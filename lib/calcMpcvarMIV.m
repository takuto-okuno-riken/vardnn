%%
% Calculate mPCVAR (multivaliate Principal Component Vector Auto-Regression) Mean Impact Value (MIV)
% returns Mean Impact Value matrix (MIV) and Mean Absolute Impact Value matrix (MAIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPCVAR network
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [MIV, MAIV] = calcMpcvarMIV(X, exSignal, nodeControl, exControl, net, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    inputNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set control 3D matrix (node x node x lags)
    [nodeControl, exControl, control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % first, calculate vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end

    MIV = nan(nodeNum,nodeMax);
    MAIV = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);

        % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        mc = net.maxComp{i};
        mu = net.mu{i};

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            Xtj1 = Yj;
            Xtj2 = Yj;
            for k=1:lags
                n = j + inputNum * (k-1);
                Xtj1(:,n) = Xtj1(:,n) * 1.1;
                Xtj2(:,n) = Xtj2(:,n) * 0.9;
            end
            scorej1 = (Xtj1(:,idx) - mu) / net.coeff{i}.';
            scorej2 = (Xtj2(:,idx) - mu) / net.coeff{i}.';
            pcXti1 = [scorej1(:,1:mc), ones(sigLen-lags,1)]; % might not be good to add bias
            pcXti2 = [scorej2(:,1:mc), ones(sigLen-lags,1)]; % might not be good to add bias
            IV1 = pcXti1 * net.bvec{i};
            IV2 = pcXti2 * net.bvec{i};

            MIV(i,j) = mean(IV1 - IV2);
            MAIV(i,j) = mean(abs(IV1 - IV2));
        end
    end
end

