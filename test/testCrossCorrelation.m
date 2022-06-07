
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
    figure; [NCC, lags] = plotCrossCorrelation(si, exSignal, [], exControl, 5); % replaced by faster version
    [NCC2, lags2] = calcCrossCorrelation(si, exSignal, [], exControl, 5); % old version
    sum(abs(NCC-NCC2),'all')

    figure; plotPartialCorrelation(si, exSignal, [], exControl);
    figure; [NCC, lags] = plotPartialCrossCorrelation(si, exSignal, [], exControl, 5);
%    [NCC2, lags] = calcPartialCrossCorrelation_(si, exSignal, [], exControl, 5); % this one is bad

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
    [NCC2, lags2] = calcCrossCorrelation(si, exSignal, [], exControl, 5, 1); % old version
    sum(abs(NCC-NCC2),'all')

    figure; plotPartialCorrelation(si, exSignal, [], exControl, 1);
    figure; [NCC, lags] = plotPartialCrossCorrelation(si, exSignal, [], exControl, 5, 1);

    %% test pattern 3
    n = 8;
    load(['results/simsc/sin-' num2str(n) 'x2000-idx6-1-1-result.mat']);
    plotTwoSignals(S(:,:,1),S(:,:,2),0,[0, 1]);
    figure; [NCC1, lags] = plotCrossCorrelation(S(:,:,1), [], [], [], 5);
    figure; [NCC2, lags] = plotCrossCorrelation(S(:,:,2), [], [], [], 5);
    [NCC3, ~] = calcCrossCorrelation(S(:,:,2), [], [], [], 5); % old version
    sum(abs(NCC3-NCC2),'all') % check diff
    U = tril(nan(n));
    NCC1(:,:,6)=[]; NCC2(:,:,6)=[];
    getCosSimilarity(NCC1+U, NCC2+U)
    figure; [PNCC1, lags] = plotPartialCrossCorrelation(S(:,:,1), [], [], [], 5);
    figure; [PNCC2, lags] = plotPartialCrossCorrelation(S(:,:,2), [], [], [], 5);
    PNCC1(:,:,6)=[]; PNCC2(:,:,6)=[];
    getCosSimilarity(PNCC1+U, PNCC2+U)
    
    load(['results/simsc/v09-' num2str(n) 'x2000-idx6-1-1-result.mat']);
    plotTwoSignals(S(:,:,1),S(:,:,2),0,[0, 1]);
    figure; [NCC1, lags] = plotCrossCorrelation(S(:,:,1), [], [], [], 5);
    figure; [NCC2, lags] = plotCrossCorrelation(S(:,:,2), [], [], [], 5);
    [NCC3, ~] = calcCrossCorrelation(S(:,:,2), [], [], [], 5); % old version
    sum(abs(NCC3-NCC2),'all') % check diff
    U = tril(nan(n));
    NCC1(:,:,6)=[]; NCC2(:,:,6)=[];
    getCosSimilarity(NCC1+U, NCC2+U)
    figure; [PNCC1, lags] = plotPartialCrossCorrelation(S(:,:,1), [], [], [], 5);
    figure; [PNCC2, lags] = plotPartialCrossCorrelation(S(:,:,2), [], [], [], 5);
    PNCC1(:,:,6)=[]; PNCC2(:,:,6)=[];
    getCosSimilarity(PNCC1+U, PNCC2+U)
end

