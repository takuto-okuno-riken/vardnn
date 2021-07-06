% Before using this function, download PartiallyConditionedGrangerCausality codes from
% https://github.com/danielemarinazzo/PartiallyConditionedGrangerCausality
% and add a path "PartiallyConditionedGrangerCausality-master" and sub folders. 

function testPCGC
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    lags = 2;
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
    figure; GC = plotPCGC(si, lags, 0.8, 0);
    
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0);
end

