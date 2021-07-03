%%
% get mVARDNN training result (mean value of each nodes)
% input:
%  net         trained VARDNN network

function [time, loss, rsme] = getMvarDnnTrainingResult(net)
    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    maxEpochs = length(net.trainInfo{1, 1}.TrainingLoss);
    a=0; b=0; c=0;
    for k=1:nodeNum
        a = a + net.trainTime;
        b = b + net.trainInfo{k, 1}.TrainingLoss(1,maxEpochs);
        c = c + net.trainInfo{k, 1}.TrainingRMSE(1,maxEpochs);
    end
    time = a / nodeNum;
    loss = b / nodeNum;
    rsme = c / nodeNum;
end