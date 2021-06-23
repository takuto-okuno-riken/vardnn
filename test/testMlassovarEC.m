
function testMlassovarEC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 3;
    lambda = 0.01;
    elaAlpha = 0.5;
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
    netMVAR = initMlassovarNetwork(si, exSignal, [], exControl, lags, lambda, elaAlpha);
    
    % show multivaliate LassoVAR-EC
    figure; pEC = plotPlassovarEC(si, exSignal, [], exControl, lags, lambda, elaAlpha, 0);
    figure; mEC = plotMlassovarEC(netMVAR, [], exControl, 0);
    % show multivaliate & pairwise LassoVAR-GC
    figure; pGC = plotPlassovarGCI(si, exSignal, [], exControl, lags, lambda, elaAlpha, 0);
    figure; mGC = plotMlassovarGCI(si, exSignal, [], exControl, netMVAR, 0);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0);

    %% test pattern 2
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % init MLassoVAR network
    netMVAR = initMlassovarNetwork(si, exSignal, [], exControl, lags, lambda, elaAlpha);
    
    % show multivaliate MVAR-EC
    figure; pEC = plotPlassovarEC(si, exSignal, [], exControl, lags, lambda, elaAlpha, 0, 1);
    figure; mEC = plotMlassovarEC(netMVAR, [], exControl, 0, 1);
    % show multivaliate & pairwise LassoVAR-GC
    figure; pGC = plotPlassovarGCI(si, exSignal, [], exControl, lags, lambda, elaAlpha, 0, 0.05, 1);
    figure; mGC = plotMlassovarGCI(si, exSignal, [], exControl, netMVAR, 0, 0.05, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0, 0, 1);
end

