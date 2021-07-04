
function testTrainDlcm
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    si = si(1:36,1:182);
%%{
    y  = si;
    width = 5;
    siend = ceil(size(si,2)/width);
    for i=1:siend
        for j=1:width
            k = (i-1)*width+j;
            y(:,k) = si(:,i)*(width-(j-1))/width + si(:,i+1)*(j-1)/width;
        end
    end
    si = y;
%%}
    % load network or training network
    netFile = ['results/dlcm-net-test' num2str(size(si,1)) 'a.mat'];
    if exist(netFile, 'file')
        load(netFile);
    else
        % init VARDNN network
        netDLCM = initMvarDnnNetwork(si);

        % set training options
        maxEpochs = 1000;
        sigLen = size(si,2);
        miniBatchSize = ceil(sigLen / 3);

        options = trainingOptions('adam', ...
            'ExecutionEnvironment','cpu', ...
            'MaxEpochs',maxEpochs, ...
            'MiniBatchSize',miniBatchSize, ...
            'Shuffle','every-epoch', ...
            'GradientThreshold',1,...
            'Verbose',false);
%            'Plots','training-progress');
%{
            'InitialLearnRate',0.001, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',400, ...
            'LearnRateDropFactor',0.5, ...
%}
        % training VARDNN network
        netDLCM = trainMvarDnnNetwork(si, [], [], [], netDLCM, options);
        save(netFile, 'netDLCM');
    end
    [time, loss, rsme] = getMvarDnnTrainingResult(netDLCM);
    disp(['train result (mean nodes) time=' num2str(time) ', loss=' num2str(loss) ', rsme=' num2str(rsme)]);
    
    % weight plot
    disp('plot connection weight amang nodes');
    nodeNum = length(netDLCM.nodeLayers);
    C = zeros(nodeNum,nodeNum);
    for i=1:nodeNum
        weight = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        mweight = mean(weight,1);
        eweight = std(weight,1) / sqrt(size(weight,1));
        C(i,:) = weight(1:nodeNum);
        % bar plot
        %{
        figure;
        bar(mweight);
        hold on
        er = errorbar(1:length(mweight),mweight,eweight);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        hold off
        %}
        % box-and-whisker plot
        % figure;
        % boxplot(weight);
    end
    % show matrix
    figure;
    image(C,'CDataMapping','scaled');
    colorbar;

    % show functional conectivity of original node signals
    figure; FC = plotFunctionalConnectivity(si);

    % test & show predicted
    figure; [S, time, mae, maeerr] = plotPredictSignals(si, [], [], [], netDLCM, 0);
    
    % show functional conectivity of predicted node signals
    figure; FC = plotFunctionalConnectivity(S);

    % plot correlation graph between original predicted node signals
    figure; R = plotTwoSignalsCorrelation(si, S) % show R result
end
