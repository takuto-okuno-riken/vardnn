%%
% Caluclate pairwise Granger Causality
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
    if nargin < 7
        isFullNode = 0;
    end
    if nargin < 6
        alpha = 0.05;
    end
    if nargin < 5
        lags = 3;
    end
    if nargin < 4
        exControl = [];
    end
    if nargin < 3
        nodeControl = [];
    end
    if nargin < 2
        exSignal = [];
    end
    nodeNum = size(X,1);
    nodeInNum = nodeNum + size(exSignal,1);
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % set node input
    if ~isempty(exSignal)
        X = [X; exSignal];
    end
    
    gcI = nan(nodeNum,nodeMax);
    h = nan(nodeNum,nodeMax);
    P = nan(nodeNum,nodeMax);
    F = nan(nodeNum,nodeMax);
    cvFd = nan(nodeNum,nodeMax);
    AIC = nan(nodeNum,nodeMax);
    BIC = nan(nodeNum,nodeMax);
    for i=1:nodeNum
        Y = X;
        if ~isempty(nodeControl)
            filter = nodeControl(i,:).';
            Y(1:nodeNum,:) = Y(1:nodeNum,:) .* filter;
        end
        if ~isempty(exControl)
            filter = exControl(i,:).';
            Y(nodeNum+1:end,:) = Y(nodeNum+1:end,:) .* filter;
        end
        for j=1:nodeMax
            if i==j, continue; end
            [gcI(i,j), h(i,j), P(i,j), F(i,j), cvFd(i,j), AIC(i,j), BIC(i,j)] = calcPairGrangerCausality(Y(i,:), Y(j,:), lags, alpha);
        end
    end
end

