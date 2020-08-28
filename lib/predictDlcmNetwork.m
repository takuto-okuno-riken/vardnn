%%
% Prot predicted signals by traind DLCM
% input:
%  X            multivariate time series matrix (node x time series)
%  inSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  netDLCM      trained DLCM network

function [S, time, mae, maeerr] = predictDlcmNetwork(X, inSignal, nodeControl, inControl, netDLCM)
    nodeNum = size(X,1);
    a=0; b=0; c=0; d=0; e=0; errs=[]; S=X;
    maxEpochs = length(netDLCM.trainInfo{1, 1}.TrainingRMSE);
    for i=1:nodeNum
        if isempty(inSignal)
            nodeInputOrg = X(:,1:end-1);
        else
            nodeInputOrg = [X(:,1:end-1); inSignal(:,1:end-1)];
        end
        nodeInput = nodeInputOrg;
        nodeTeach = single(X(i,2:end));
        if ~isempty(nodeControl)
            filter = repmat(nodeControl(i,:).', 1, size(nodeControl,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        if ~isempty(inControl)
            filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        zPred = predict(netDLCM.nodeNetwork{i}, nodeInput);
        S(i,2:end) = zPred;

        d = d + sqrt(immse(zPred, nodeTeach));
        e = e + mean(abs(zPred - nodeTeach));
        a = a + netDLCM.trainTime;
        b = b + netDLCM.trainInfo{i, 1}.TrainingLoss(1,maxEpochs);
        c = c + netDLCM.trainInfo{i, 1}.TrainingRMSE(1,maxEpochs);
        A = single(zPred - nodeTeach);
        errs = [errs, A];
    end
    time = a / nodeNum;
    mae = e / nodeNum;
    maeerr = std(errs) / sqrt(length(errs));
end
