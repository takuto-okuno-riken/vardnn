
function performanceCheckSignalLen
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    
    nodeNum = 32;
    inputNum = 10;
    hiddenNums = [48; ceil(48*2/3)];

    start = 100;
    step = 100;
    stepMax = 20;
    for i=start:step:step*stepMax
        si = siOrg(1:nodeNum,1:i);
        inSignal = [];
        if inputNum > 0
            inSignal = zeros(inputNum, size(si,2)); % exogenous input matrix
        end
        sigLen = size(si,2);

        % set training options
        maxEpochs = 500;
        miniBatchSize = ceil(sigLen / 3);

        options = trainingOptions('adam', ...
            'ExecutionEnvironment','cpu', ...
            'MaxEpochs',maxEpochs, ...
            'MiniBatchSize',miniBatchSize, ...
            'Shuffle','every-epoch', ...
            'GradientThreshold',1,...
            'Verbose',false);
    %            'Plots','training-progress');

        disp(['training signal length ' num2str(i)]);
        % performance check of hidden layers
        dlcmFile = ['results/net-siglen-' num2str(i) '-hdn' num2str(hiddenNums(1)) '-' num2str(hiddenNums(2)) '.mat'];
        if exist(dlcmFile, 'file')
            continue;
        end

        % layer parameters
        netDLCM = createDlcmNetwork(nodeNum, inputNum, hiddenNums);

        % training DLCM network
        netDLCM = trainDlcmNetwork(si, inSignal, [], [], netDLCM, options);
        save(dlcmFile, 'netDLCM');
    end
    
    % get result
    resultTime = zeros(stepMax,1);
    resultLoss = zeros(stepMax,1);
    resultRSME = zeros(stepMax,1);
    resultRSMEAll = zeros(stepMax,1);
    resultMAEAll = zeros(stepMax,1);
    resultMAEErr = zeros(stepMax,1);
    resultErrs = zeros(stepMax,nodeNum,step*stepMax);

    disp('checking prediction.');
    for i=1:stepMax
        disp(['loading result of ' num2str(i*step)]);
        % performance check of signal length
        dlcmFile = ['results/net-siglen-' num2str(i*step) '-hdn' num2str(hiddenNums(1)) '-' num2str(hiddenNums(2)) '.mat'];
        load(dlcmFile);

        % set signals
        si = siOrg(1:nodeNum,1:i*step);
        inSignal = [];
        if inputNum > 0
            inSignal = zeros(inputNum, size(si,2)); % exogenous input matrix
        end
        
        a=0; b=0; c=0; d=0; e=0;
        for k=1:nodeNum
            if isempty(inSignal)
                nodeInput = si(:,1:end-1);
            else
                nodeInput = [si(:,1:end-1); inSignal(:,1:end-1)];
            end
            nodeTeach = single(si(k,2:end));
            zPred = predict(netDLCM.nodeNetwork{k}, nodeInput);
            d = d + sqrt(immse(zPred, nodeTeach));
            e = e + mean(abs(zPred - nodeTeach));
            a = a + netDLCM.trainTime;
            b = b + netDLCM.trainInfo{k, 1}.TrainingLoss(1,maxEpochs);
            c = c + netDLCM.trainInfo{k, 1}.TrainingRMSE(1,maxEpochs);
            resultErrs(i,k,1:i*step-1) = single(zPred - nodeTeach);
        end

        resultTime(i,1) = a / nodeNum;
        resultLoss(i,1) = b / nodeNum;
        resultRSME(i,1) = c / nodeNum;
        resultRSMEAll(i,1) = d / nodeNum;
        resultMAEAll(i,1) = e / nodeNum;
        A = resultErrs(i,:,1:i*step-1);
        errs = A(:);
        resultMAEErr(i,1) = std(errs,1) / sqrt(length(errs)); % uncorrected std error
    end
    
    % show error bar graph
    figure;
    errorbar(start:step:step*stepMax, resultMAEAll, resultMAEErr);

    % seve result file
    filename = 'results/result-signal-length.mat';
    save(filename, 'resultTime','resultLoss','resultRSME','resultRSMEAll','resultMAEAll','resultMAEErr','resultErrs');
end
