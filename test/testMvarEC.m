
function testMvarEC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 5;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1);

    %% test pattern 1 
    netMVAR = initMvarNetwork(si, exSignal, [], exControl, lags);
    
    % show multivaliate MVAR-EC
    figure; mEC = plotMvarEC(netMVAR, [], exControl, 0);
    figure; pEC = plotPvarEC(si, exSignal, [], exControl, lags, 0);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0);

    %% test pattern 2
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % init MVAR network
    netMVAR = initMvarNetwork(si, exSignal, [], exControl, lags);
    
    % show multivaliate MVAR-EC
    figure; mEC = plotMvarEC(netMVAR, [], exControl, 0, 1);
    figure; pEC = plotPvarEC(si, exSignal, [], exControl, lags, 0, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0, 0, 1);
end

