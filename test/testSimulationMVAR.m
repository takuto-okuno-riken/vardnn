
function testSimulationMVAR
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 4;
    lags = 5;
    sigLen = 100;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = logical(ones(nodeNum,exNum));
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);

    %% test pattern 1 
    netMVAR = initMvarNetwork(si, exSignal, [], exControl, lags);
    
    % simulate mVAR network with 1st frame & exogenous input signal
    [S, time] = simulateMvarNetwork(si, exSignal, [], exControl, netMVAR);

    figure; [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);
    % show multivaliate MVAR-DI
    figure; DI = plotMvarDI(netMVAR, [], exControl, 0);
    
    %% test pattern along with testSimulationMpcvarDnn.m
    % load signals
    load('data/marmoset-aneth-sample2-roi225.mat');
    siOrg = convert2SigmoidSignal(si);
    load('test/testTrain-rand500-uniform.mat');
    uuOrg = si;
    
    nodeNum = 100;
    exNum = 0;
    sigLen = 50;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    exControl = [];
    
    %% simulate test without exogenous, change lags
    for lags=1:3
        [time, mae, errs] = trainAndPlotMvar(si, exSignal, [], exControl, lags);
    end
        
    %% simulate test with exogenous, change lags
    exNum = nodeNum;
    exSignal = uuOrg(1:exNum,1:sigLen);
    exControl = logical(eye(nodeNum,nodeNum));
    for lags=1:3
        [time, mae, errs] = trainAndPlotMvar(si, exSignal, [], exControl, lags);
    end
    
    %% simulate test with exogenous, change lags
    exNum = 4;
    exSignal = uuOrg(1:exNum,1:sigLen);
    exControl = logical(ones(nodeNum,exNum));
    for lags=1:3
        [time, mae, errs] = trainAndPlotMvar(si, exSignal, [], exControl, lags);
    end
    
    %% simulate test with exogenous, change node num and sigLen
    nodeNums = [100,200,300,400];
    sigLens = [50,100,150,200,250,300,350];
    allErrs = nan(max(sigLens), 5*length(nodeNums)*length(sigLens));
    for i=1:length(nodeNums)
        nodeNum = nodeNums(i);
        for j=1:length(sigLens)
            sigLen = sigLens(j);
            for lags=1:5
                si = siOrg(1:nodeNum,1:sigLen);
                exSignal = uuOrg(1:exNum,1:sigLen);
                exControl = logical(ones(nodeNum,exNum));
                [time, mae, errs] = trainAndPlotMvar(si, exSignal, [], exControl, lags);
                % set to one matrix
                if size(errs,3) > 1
                    errX = squeeze(mean(abs(errs),1));
                else
                    errX = mean(abs(errs),1).';
                end
                errXY = mean(errX,1);
                [M, I] = min(errXY);
                allErrs(1:sigLen, (j-1)*(5*4)+(lags-1)*4+i) = errX(:,I);
            end
        end
    end
end


function [time, mae, errs] = trainAndPlotMvar(si, exSignal, nodeControl, exControl, lags)
    nodeNum = size(si,1);
    exNum = size(exSignal,1);
    sigLen = size(si,2);
    resFile = ['results/mvar/sim-test-' num2str(nodeNum) '_' num2str(exNum) 'x' num2str(sigLen) '-l' num2str(lags) '.mat'];
    if exist(resFile, 'file')
        load(resFile);
    else
        net = initMvarNetwork(si, exSignal, nodeControl, exControl, lags);
        % simulate mVAR network with 1st frame & exogenous input signal
        [S, time] = simulateMvarNetwork(si, exSignal, nodeControl, exControl, net);
        [mae, ~, errs] = getTwoSignalsError(si, S);
        % recover training 
        [trainedNet, time, mae, errs] = recoveryTrainMvarNetwork(si, exSignal, nodeControl, exControl, net, 0.008, [0.1, 0.4]);
        save(resFile, 'time', 'mae', 'errs', 'S');
    end
    disp(['train result time=' num2str(time) ', mae=' num2str(mae)]);
    errX = squeeze(mean(abs(errs),1));
    figure; plot(errX);
    title([num2str(nodeNum) '_' num2str(exNum) 'x' num2str(sigLen) ' lags=' num2str(lags)]); drawnow;
end

