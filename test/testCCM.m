% Before using this function, download xmap codes from
% https://github.com/danm0nster/xmap
% and add a path "xmap-master" folder.

function testCCM
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

    %% test pattern 1 
    for lags=1:3
        % show Convergent Cross Mapping
        figure; CCM = plotConvCrossMap(si, exSignal, [], exControl, lags);
        CCM2 = calcConvCrossMap(si, exSignal, [], exControl, lags);
        if ~isequaln(CCM,CCM2)
            disp('error : CCM1 != CCM2 !');
            return;
        end
        % compare to mvGC
        figure; GC = plotMultivariateGCI(si, exSignal, [], exControl, lags, 0);
        figure; GC = plotPairwiseGCI(si, exSignal, [], exControl, lags, 0);
    end

    %% test pattern 2
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    % show Convergent Cross Mapping
    figure; CCM = plotConvCrossMap(si, exSignal, [], exControl, lags, 1, 1);
    CCM2 = calcConvCrossMap(si, exSignal, [], exControl, lags, 1, [], 'linear', 1);
    if ~isequaln(CCM,CCM2)
        disp('error : CCM1 != CCM2 !');
        return;
    end
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, nodeControl, exControl, lags, 0, 0, 1);
    figure; GC = plotPairwiseGCI(si, exSignal, nodeControl, exControl, lags, 0, 0, 1);

    %% test pattern 3
    lags = 3;
    exNum = 4;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    nodeControl = ones(nodeNum,nodeNum,lags);
    for k=1:nodeNum, nodeControl(k,k,2)=0; end
    exControl = ones(nodeNum,exNum,lags);
    si(3,2:end) = exSignal(1,1:sigLen-1);
    si(5,3:end) = si(1,1:sigLen-2); % lag=2, this will be blocked by nodeControl
    nodeControl(5,1,2) = 0; % <= comment out and check control effect
    si(7,3:end) = exSignal(2,1:sigLen-2); % lag=2, this will be blocked by exControl
    exControl(7,2,2) = 0; % <= comment out and check control effect

    % show Convergent Cross Mapping (3D control does not work)
    figure; CCM = plotConvCrossMap(si, exSignal, nodeControl, exControl, lags, 1, 1);
    % compare to mvGC
    figure; GC = plotMultivariateGCI(si, exSignal, nodeControl, exControl, lags, 0, 0, 1);
    figure; GC = plotPairwiseGCI(si, exSignal, nodeControl, exControl, lags, 0, 0, 1);
end

