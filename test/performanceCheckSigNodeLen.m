
function performanceCheckSignalLen
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    
    exNum = 0;
    hiddenNums = [48; ceil(48*2/3)];

    start = 100;
    step = 100;
    stepMax = 12;
    start2 = 15;
    step2 = 15;
    stepMax2 = 12;
    for i=start:step:step*stepMax
        for j=start2:step2:step2*stepMax2
            si = siOrg(1:j,1:i);
            nodeNum = size(si,1);
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

            disp(['training signal ' num2str(j) 'x' num2str(i)]);
            % performance check of hidden layers
            dlcmFile = ['results/net-sig' num2str(j) 'x' num2str(i) '-hdn' num2str(hiddenNums(1)) '-' num2str(hiddenNums(2)) '.mat'];
            if exist(dlcmFile, 'file')
                continue;
            end

            % layer parameters
            netDLCM = createDlcmNetwork(nodeNum, exNum, hiddenNums);

            % training DLCM network
            netDLCM = trainDlcmNetwork(si, [], [], [], netDLCM, options);
            save(dlcmFile, 'netDLCM');
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
            disp(['loading result of ' num2str(j*step2) 'x' num2str(i*step)]);
            % performance check of signal length
            dlcmFile = ['results/net-sig' num2str(j*step2) 'x' num2str(i*step) '-hdn' num2str(hiddenNums(1)) '-' num2str(hiddenNums(2)) '.mat'];
            load(dlcmFile);

            % set signals
            si = siOrg(1:j*step2,1:i*step);
            exSignal = [];
            if exNum > 0
                exSignal = zeros(exNum, size(si,2)); % exogenous input matrix
            end
            nodeNum = size(si,1);
            sigLen = size(si,2);
        
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
    %image(resultRSMEAll,'CDataMapping','scaled');
    colorbar;

    % show matrix
    figure;
    imagesc(resultMAEAll,clims);
    colorbar;

    % seve result file
    filename = 'results/result-sig-nodelen.mat';
    save(filename, 'resultTime','resultLoss','resultRSME','resultRSMEAll','resultMAEAll','resultMAEErr');
end
