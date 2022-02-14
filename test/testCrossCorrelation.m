
function testCrossCorrelation
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 1;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    nodeControl = [];
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1); % caution! node 2 & 4 is Multicollinearity case (correlated)
    si(1,2:end) = si(4,1:sigLen-1);

    %% test pattern 1 
    figure; plotFunctionalConnectivity(si, exSignal, [], exControl);
    figure; plotPairwiseGCI(si, exSignal, [], exControl, 3);
    figure; [NCC, lags] = plotCrossCorrelation(si, exSignal, [], exControl, 5);

    figure; plotPartialCorrelation(si, exSignal, [], exControl);
    figure; [NCC, lags] = plotPartialCrossCorrelation(si, exSignal, [], exControl, 5);
    
    %% test pattern 2
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    figure; FC = plotFunctionalConnectivity(si, exSignal, [], exControl, 1);
    figure; pGC = plotPairwiseGCI(si, exSignal, [], exControl, lags, 10, 0.05, 1);
    figure; [NCC, lags] = plotCrossCorrelation(si, exSignal, [], exControl, 5, 1);

    figure; plotPartialCorrelation(si, exSignal, [], exControl, 1);
    figure; [NCC, lags] = plotPartialCrossCorrelation(si, exSignal, [], exControl, 5, 1);
end

