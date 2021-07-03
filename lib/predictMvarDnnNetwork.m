%%
% Prot predicted signals by traind mVAR DNN
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained mVAR DNN network

function [Y, time] = predictMvarDnnNetwork(X, exSignal, nodeControl, exControl, net)
    nodeNum = size(X,1);
    Y = X;
    if isempty(exSignal)
        nodeInputOrg = X;
    else
        nodeInputOrg = [X; exSignal];
    end
    ticH = tic;
    for i=1:nodeNum
        nodeInput = nodeInputOrg;
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', 1, size(nodeInput,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        if ~isempty(exControl)
            filter = repmat(exControl(i,:).', 1, size(nodeInput,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        zPred = predict(net.nodeNetwork{i}, nodeInput);
        Y(i,:) = zPred;
    end
    time = toc(ticH);
end
