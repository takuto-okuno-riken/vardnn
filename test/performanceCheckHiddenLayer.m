
function performanceCheckHiddenLayer
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    si = si(1:32,1:200);

    exNum = 10;
    exSignal = [];
    if exNum > 0
        exSignal = zeros(exNum, size(si,2)); % exogenous input matrix
    end
    
    sigLen = size(si,2);
    nodeNum = size(si,1);

    % set training options
    maxEpochs = 500;
    miniBatchSize = ceil(sigLen / 3);

    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',1,...
        'Verbose',false);
%            'Plots','training-progress');

    start = 4;
    step = 4;
    stepMax = 12;
    for i=start:step:step*stepMax
        for j=start:step:step*stepMax
            disp(['training hidden layer num ' num2str(i) '-' num2str(j)]);
            % performance check of hidden layers
            netFile = ['results/net-hidden-' num2str(i) '-' num2str(j) '.mat'];
            if exist(netFile, 'file')
                continue;
            end

            % layer parameters
            hiddenNums = [i; j];
            netDLCM = createMvarDnnNetwork(nodeNum, exNum, hiddenNums);

            % training VARDNN network
            netDLCM = trainMvarDnnNetwork(si, exSignal, [], [], netDLCM, options);
            save(netFile, 'netDLCM');
        end
    end
    
    % get result
    resultTime = zeros(stepMax,stepMax);
    resultLoss = zeros(stepMax,stepMax);
    resultRSME = zeros(stepMax,stepMax);
    resultRSMEAll = zeros(stepMax,stepMax);
    resultMAEAll = zeros(stepMax,stepMax);

    disp('checking prediction.');
    for i=1:stepMax
        for j=1:stepMax
            disp(['loading result of ' num2str(i*step) '-' num2str(j*step)]);
            % performance check of hidden layers
            netFile = ['results/net-hidden-' num2str(i*step) '-' num2str(j*step) '.mat'];
            load(netFile);

            a=0; b=0; c=0; d=0; e=0;
            for k=1:nodeNum
                if isempty(exSignal)
                    nodeInput = si(:,1:end-1);
                else
                    nodeInput = [si(:,1:end-1); exSignal(:,1:end-1)];
                end
                nodeTeach = single(si(k,2:end));
                zPred = predict(netDLCM.nodeNetwork{k}, nodeInput);
                d = d + sqrt(immse(zPred, nodeTeach));
                e = e + mean(abs(zPred - nodeTeach));
                a = a + netDLCM.trainTime;
                b = b + netDLCM.trainInfo{k, 1}.TrainingLoss(1,maxEpochs);
                c = c + netDLCM.trainInfo{k, 1}.TrainingRMSE(1,maxEpochs);
            end

            resultTime(i,j) = a / nodeNum;
            resultLoss(i,j) = b / nodeNum;
            resultRSME(i,j) = c / nodeNum;
            resultRSMEAll(i,j) = d / nodeNum;
            resultMAEAll(i,j) = e / nodeNum;
        end
    end
    % show matrix
    figure;
    image(resultLoss,'CDataMapping','scaled');
    colorbar;

    % show matrix
    figure;
    image(resultRSMEAll,'CDataMapping','scaled');
    colorbar;

    % show matrix
    figure;
    image(resultMAEAll,'CDataMapping','scaled');
    colorbar;

    % seve result file
    filename = 'results/result-hidden-check.mat';
    save(filename, 'resultTime','resultLoss','resultRSME','resultRSMEAll','resultMAEAll');
end
