
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
    % show multivaliate MVAR-EC
    figure; EC = plotMvarEC(netMVAR, [], exControl, 0);
end

