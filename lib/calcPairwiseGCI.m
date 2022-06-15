%%
% Calculate pairwise Granger Causality
% returns Granger causality index (gcI), significance (h=1 or 0)
% p-values (P), F-statistic (F) and the critical value from the F-distribution (cvFd)
% and AIC, BIC
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for autoregression (default:3)
%  alpha        the significance level of F-statistic (default:0.05)
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [gcI, h, P, F, cvFd, AIC, BIC] = calcPairwiseGCI(X, exSignal, nodeControl, exControl, lags, alpha, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, alpha = 0.05; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);
    
    gcI = nan(nodeNum,nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            [gcI(i,j), h(i,j), P(i,j), F(i,j), cvFd(i,j), AIC(i,j), BIC(i,j)] = calcPairGrangerCausality(Y(i,:), Y(j,:), lags, alpha, control(i,i,:),control(i,j,:));
        end
    end
    % output control
    if isFullNode==0
        gcI = gcI(:,1:nodeNum);
        F = F(:,1:nodeNum);
        P = P(:,1:nodeNum);
        cvFd = cvFd(:,1:nodeNum);
        h = h(:,1:nodeNum);
        AIC = AIC(:,1:nodeNum);
        BIC = BIC(:,1:nodeNum);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl(:,:,1)); nodeControl(nodeControl==0) = nan;
        gcI(:,1:nodeNum) = gcI(:,1:nodeNum) .* nodeControl;
        F(:,1:nodeNum) = F(:,1:nodeNum) .* nodeControl;
        P(:,1:nodeNum) = P(:,1:nodeNum) .* nodeControl;
        cvFd(:,1:nodeNum) = cvFd(:,1:nodeNum) .* nodeControl;
        h(:,1:nodeNum) = h(:,1:nodeNum) .* nodeControl;
        AIC(:,1:nodeNum) = AIC(:,1:nodeNum) .* nodeControl;
        BIC(:,1:nodeNum) = BIC(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl(:,:,1)); exControl(exControl==0) = nan;
        gcI(:,nodeNum+1:end) = gcI(:,nodeNum+1:end) .* exControl;
        F(:,nodeNum+1:end) = F(:,nodeNum+1:end) .* exControl;
        P(:,nodeNum+1:end) = P(:,nodeNum+1:end) .* exControl;
        cvFd(:,nodeNum+1:end) = cvFd(:,nodeNum+1:end) .* exControl;
        h(:,nodeNum+1:end) = h(:,nodeNum+1:end) .* exControl;
        AIC(:,nodeNum+1:end) = AIC(:,nodeNum+1:end) .* exControl;
        BIC(:,nodeNum+1:end) = BIC(:,nodeNum+1:end) .* exControl;
    end
end

