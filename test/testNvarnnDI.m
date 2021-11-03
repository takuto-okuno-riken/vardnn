
function testNvarnnDI
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    lags = 1;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);

    % set training options
    maxEpochs = 2000;
    sigLen = size(si,2);
    miniBatchSize = ceil(sigLen / 3);

    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'Shuffle','every-epoch', ...
        'L2Regularization',0.01, ...
        'GradientThreshold',5,...
        'Verbose',true);
%            'Plots','training-progress');

    %% test pattern 1 
    % do training or load nVARNN network
    for i=1:5
        netFile = ['results/nvarnn' num2str(i) '-di-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            % init nVARNN network
            net = initNvarnnNetwork(si, exSignal, [], exControl, i);
            % training nVARNN network
            net = trainNvarnnNetwork(si, exSignal, [], exControl, net, options);

            loss = net.trainInfo.TrainingLoss(1,maxEpochs);
            rsme = net.trainInfo.TrainingRMSE(1,maxEpochs);
            disp(['train result time=' num2str(net.trainTime) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(netFile, 'net');
        end

        % show nVARNN-GC
        figure; DI = plotNvarnnDI(net, [], exControl, 0);
        figure; MIV = plotNvarnnMIV(si, exSignal, [], exControl, net, 0);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], [], i, 0);
    end

    %% test pattern 2 -- exogenous signals
    % do training or load nVARNN network
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    nodeControl = ones(1,nodeNum);
    exControl = ones(1,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    for i=1:5
        netFile = ['results/nvarnn' num2str(i) '-di2-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            % init nVARNN network
            net = initNvarnnNetwork(si, exSignal, nodeControl, exControl, i);
            % training nVARNN network
            net = trainNvarnnNetwork(si, exSignal, nodeControl, exControl, net, options);

            loss = net.trainInfo.TrainingLoss(1,maxEpochs);
            rsme = net.trainInfo.TrainingRMSE(1,maxEpochs);
            disp(['train result time=' num2str(net.trainTime) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(netFile, 'net');
        end

        % show nVARNN-GC
        figure; DI = plotNvarnnDI(net, nodeControl, exControl, 0, 1);
        figure; MIV = plotNvarnnMIV(si, exSignal, nodeControl, exControl, net, 0, 1);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], [], i, 0, 0, 1);
    end
    
    %% test pattern 3 -- no activation function
    % do training or load nVARNN network
    exNum = 3;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    nodeControl = ones(1,nodeNum);
    exControl = ones(1,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    for i=1:5
        netFile = ['results/nvarnn' num2str(i) '-di3-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            % init nVARNN network
            net = initNvarnnNetwork(si, exSignal, nodeControl, exControl, i);
            % training nVARNN network
            net = trainNvarnnNetwork(si, exSignal, nodeControl, exControl, net, options);

            loss = net.trainInfo.TrainingLoss(1,maxEpochs);
            rsme = net.trainInfo.TrainingRMSE(1,maxEpochs);
            disp(['train result time=' num2str(net.trainTime) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(netFile, 'net');
        end

        % show nVARNN-GC
        figure; DI = plotNvarnnDI(net, nodeControl, exControl, 0, 1);
        figure; MIV = plotNvarnnMIV(si, exSignal, nodeControl, exControl, net, 0, 1);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], [], i, 0, 0, 1);
    end

    %% test pattern 4 -- exogenous signals
    % do training or load nVARNN network
    exNum = 4;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    nodeControl = [];
    exControl = [];
    si(3,2:end) = exSignal(1,1:sigLen-1);
    for i=1:5
        if i>1
            nodeControl = ones(1,nodeNum,i);
            exControl = ones(1,exNum,i);
            si(5,3:end) = si(1,1:sigLen-2); % lag=2, this will be blocked by nodeControl
%            nodeControl(1,1,2) = 0; % <= comment out and check control effect
            si(7,3:end) = exSignal(2,1:sigLen-2); % lag=2, this will be blocked by exControl
%            exControl(1,2,2) = 0; % <= comment out and check control effect
        end
        netFile = ['results/nvarnn' num2str(i) '-di4-test' num2str(nodeNum) '-' num2str(exNum) '.mat'];
        if exist(netFile, 'file')
            load(netFile);
        else
            % init nVARNN network
            net = initNvarnnNetwork(si, exSignal, nodeControl, exControl, i);
            % training nVARNN network
            net = trainNvarnnNetwork(si, exSignal, nodeControl, exControl, net, options);

            loss = net.trainInfo.TrainingLoss(1,maxEpochs);
            rsme = net.trainInfo.TrainingRMSE(1,maxEpochs);
            disp(['train result time=' num2str(net.trainTime) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
            save(netFile, 'net');
        end

        % show nVARNN-GC
        figure; DI = plotNvarnnDI(net, nodeControl, exControl, 0, 1);
        figure; MIV = plotNvarnnMIV(si, exSignal, nodeControl, exControl, net, 0, 1);
        figure; gcI = plotMultivariateGCI(si, exSignal, [], [], i, 0, 0, 1);
    end
end

