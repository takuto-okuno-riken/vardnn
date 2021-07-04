
function performanceCheckVariHidden
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    
    exNum = 0;

    start = 200;
    step = 200;
    stepMax = 4;
    start2 = 100;
    step2 = 100;
    stepMax2 = 4;
    for i=start:step:step*stepMax
        for j=start2:step2:step2*stepMax2
            si = siOrg(1:j,1:i);
            nodeNum = size(si,1);
            sigLen = size(si,2);

            hdnNums = estimateHiddenNeurons(nodeNum, sigLen);

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

            disp(['training signal ' num2str(j) 'x' num2str(i)]);
            % performance check of hidden layers
            netFile = ['results/net-vari' num2str(j) 'x' num2str(i) '-hdn' num2str(hdnNums(1)) '-' num2str(hdnNums(2)) '.mat'];
            if exist(netFile, 'file')
                continue;
            end

            % layer parameters
            netDLCM = createMvarDnnNetwork(nodeNum, exNum, hdnNums);

            % training VARDNN network
            netDLCM = trainMvarDnnNetwork(si, [], [], [], netDLCM, options);
            save(netFile, 'netDLCM');
        end
    end
    
    % get result
    resultTime = zeros(stepMax,stepMax2);
    resultLoss = zeros(stepMax,stepMax2);
    resultRSME = zeros(stepMax,stepMax2);
    resultRSMEAll = zeros(stepMax,stepMax2);
    resultMAEAll = zeros(stepMax,stepMax2);
    resultMAEErr = zeros(stepMax,stepMax2);

    disp('checking prediction.');
    for i=1:stepMax
        for j=1:stepMax2
            nodeNum = j*step2;
            sigLen = i*step;
            hdnNums = estimateHiddenNeurons(nodeNum, sigLen);

            disp(['loading result of ' num2str(nodeNum) 'x' num2str(sigLen)]);
            % performance check of signal length
            netFile = ['results/net-vari' num2str(nodeNum) 'x' num2str(sigLen) '-hdn' num2str(hdnNums(1)) '-' num2str(hdnNums(2)) '.mat'];
            load(netFile);

            % set signals
            si = siOrg(1:nodeNum,1:sigLen);
            exSignal = [];
            if exNum > 0
                exSignal = zeros(exNum, size(si,2)); % exogenous input matrix
            end
        
            a=0; b=0; c=0; d=0; e=0; errs=[];
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
                A = single(zPred - nodeTeach);
                errs = [errs, A];
            end

            resultTime(j,i) = a / nodeNum;
            resultLoss(j,i) = b / nodeNum;
            resultRSME(j,i) = c / nodeNum;
            resultRSMEAll(j,i) = d / nodeNum;
            resultMAEAll(j,i) = e / nodeNum;
            resultMAEErr(j,i) = std(errs,1) / sqrt(length(errs)); % uncorrected std error
        end
    end
    
    % show matrix
    figure;
    clims = [0.005 0.2];
    imagesc(resultRSMEAll,clims);
    colorbar;

    % show matrix
    figure;
    imagesc(resultMAEAll,clims);
    %image(resultMAEAll,'CDataMapping','scaled');
    colorbar;

    % seve result file
    filename = 'results/result-vari-hidden.mat';
    save(filename, 'resultTime','resultLoss','resultRSME','resultRSMEAll','resultMAEAll','resultMAEErr');
end
