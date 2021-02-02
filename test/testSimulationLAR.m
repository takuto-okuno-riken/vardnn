
function testSimulationLAR
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
    netLAR = initLarNetwork(si, exSignal, [], exControl, lags);
    
    % simulate LAR network with 1st frame & exogenous input signal
    [S, time] = simulateLarNetwork(si, exSignal, [], exControl, netLAR);

    figure; [mae, maeerr] = plotTwoSignals(si, S);
    disp(['simulation time=' num2str(time) ', mae=' num2str(mae)]);
    
    % show original & simulated signal FC
    figure; FC = plotFunctionalConnectivity(si);
    figure; FC = plotFunctionalConnectivity(S);
    % show original & simulated signal granger causality index (gc-EC)
    figure; gcI = plotPairwiseGCI(si);
    figure; gcI = plotPairwiseGCI(S);
    % show multivaliate LAR-EC
    figure; EC = plotLarEC(netLAR, [], exControl, 0);
end

