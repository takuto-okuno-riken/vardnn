%%
% get mVARDNN training result (mean value of each nodes)
% input:
%  net         trained VARDNN network

function [time, loss, rsme] = getMvarDnnTrainingResult(net)
    if isfield(net, 'nodeNum'), nodeNum = net.nodeNum; else nodeNum = length(net.nodeNetwork); end % for compatibility
    a=0; b=0; c=0; count=0;
    for k=1:nodeNum
        if isempty(net.trainInfo{k, 1}), continue; end
        maxEpochs = length(net.trainInfo{k, 1}.TrainingLoss);
        a = a + net.trainTime;
        b = b + net.trainInfo{k, 1}.TrainingLoss(1,maxEpochs);
        c = c + net.trainInfo{k, 1}.TrainingRMSE(1,maxEpochs);
        count = count + 1;
    end
    time = a / count;
    loss = b / count;
    rsme = c / count;
end