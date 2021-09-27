
function testMpcvarDI
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
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
    lags = 1;
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
    lags = 3;
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

    %% test pattern 3
    exNum = 4;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    nodeControl = ones(nodeNum,nodeNum,lags);
    for i=1:nodeNum, nodeControl(i,i,2)=0; end
    exControl = ones(nodeNum,exNum,lags);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    si(5,3:end) = si(1,1:sigLen-2); % lag=2, this will be blocked by nodeControl
    nodeControl(5,1,2) = 0; % <= comment out and check control effect
    si(7,3:end) = exSignal(2,1:sigLen-2); % lag=2, this will be blocked by exControl
    exControl(7,2,2) = 0; % <= comment out and check control effect

    % init PCVAR network
    net = initMpcvarNetwork(si, exSignal, nodeControl, exControl, lags);
    net2 = initPpcvarNetwork(si, exSignal, nodeControl, exControl, lags);
    
    % show multivaliate & pairwise MVAR-DI
    figure; mDI = plotMpcvarDI(net, nodeControl, exControl, 0, 1);
    figure; pDI = plotPpcvarDI(net2, nodeControl, exControl, 0, 1);
    % show multivaliate & pairwise PCVAR-GC
    figure; mGC = plotMpcvarGCI(si, exSignal, nodeControl, exControl, net, 0, 0.05, 1);
    figure; pGC = plotPpcvarGCI(si, exSignal, nodeControl, exControl, net2, 0, 0.05, 1);
    figure; mMIV = plotMpcvarMIV(si, exSignal, nodeControl, exControl, net, 0, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, nodeControl, exControl, lags, 0, 0.05, 1);
end

