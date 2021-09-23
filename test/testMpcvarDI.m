
function testMpcvarDI
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 3;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1); % caution! node 2 & 4 is Multicollinearity case (correlated)

    %% test pattern 1 
    net = initMpcvarNetwork(si, exSignal, [], exControl, lags);
    net2 = initPpcvarNetwork(si, exSignal, [], exControl, lags);
    
    % show multivaliate & pairwise PCVAR-DI
    figure; mDI = plotMpcvarDI(net, [], exControl, 0);
    figure; pDI = plotPpcvarDI(net2, [], exControl, 0);
    % show multivaliate & pairwise PCVAR-GC
    figure; mGC = plotMpcvarGCI(si, exSignal, [], exControl, net, 0);
    figure; pGC = plotPpcvarGCI(si, exSignal, [], exControl, net2, 0);
    figure; mMIV = plotMpcvarMIV(si, exSignal, [], exControl, net, 0);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0);

    %% test pattern 2
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % init PCVAR network
    net = initMpcvarNetwork(si, exSignal, [], exControl, lags);
    net2 = initPpcvarNetwork(si, exSignal, [], exControl, lags);
    
    % show multivaliate & pairwise MVAR-DI
    figure; mDI = plotMpcvarDI(net, [], exControl, 0, 1);
    figure; pDI = plotPpcvarDI(net2, [], exControl, 0, 1);
    % show multivaliate & pairwise PCVAR-GC
    figure; mGC = plotMpcvarGCI(si, exSignal, [], exControl, net, 0, 0.05, 1);
    figure; pGC = plotPpcvarGCI(si, exSignal, [], exControl, net2, 0, 0.05, 1);
    figure; mMIV = plotMpcvarMIV(si, exSignal, [], exControl, net, 0, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0, 0.05, 1);
end

