%%
% get DLCM training result (mean value of each nodes)
% input:
%  netDLCM     trained DLCM network

function [time, loss, rsme] = getDlcmTrainingResult(netDLCM)
    nodeNum = length(netDLCM.nodeNetwork);
    maxEpochs = length(netDLCM.trainInfo{1, 1}.TrainingLoss);
    a=0; b=0; c=0;
    for k=1:nodeNum
        a = a + netDLCM.trainTime;
        b = b + netDLCM.trainInfo{k, 1}.TrainingLoss(1,maxEpochs);
        c = c + netDLCM.trainInfo{k, 1}.TrainingRMSE(1,maxEpochs);
    end
    time = a / nodeNum;
    loss = b / nodeNum;
    rsme = c / nodeNum;
end