
function testSvCrossCorrelation
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
    figure; plotPartialCorrelation(si, exSignal, [], exControl);
    figure; [NCC, lags] = plotPartialCrossCorrelation(si, exSignal, [], exControl, 5);
    figure; [NCC, lags] = plotSvPartialCrossCorrelation(si, exSignal, [], exControl, 5);
    figure; [NCC, lags] = plotSvPartialCrossCorrelation(si, exSignal, [], exControl, 5, 'gaussian');
    
    %% test pattern 2
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum);
    si(3,2:end) = exSignal(1,1:sigLen-1);

    figure; plotPartialCorrelation(si, exSignal, [], exControl, 1);
    figure; [NCC, lags] = plotPartialCrossCorrelation(si, exSignal, [], exControl, 5, 1);
    figure; [NCC, lags] = plotSvPartialCrossCorrelation(si, exSignal, [], exControl, 5, 'linear', 'auto', 1);
    figure; [NCC, lags] = plotSvPartialCrossCorrelation(si, exSignal, [], exControl, 5, 'gaussian', 'auto', 1);

    %% test pattern 3
    n = 8;
    U = tril(nan(n));
    load(['results/simsc/sin-' num2str(n) 'x2000-idx6-1-1-result.mat']);
    plotTwoSignals(S(:,:,1),S(:,:,2),0,[0, 1]);
    figure; [NCC1, lags] = plotPartialCrossCorrelation(S(:,:,1), [], [], [], 5);
    figure; [NCC2, lags] = plotPartialCrossCorrelation(S(:,:,2), [], [], [], 5);
    NCC1(:,:,6)=[]; NCC2(:,:,6)=[];
    getCosSimilarity(NCC1+U, NCC2+U)
    figure; [PNCC1, lags] = plotSvPartialCrossCorrelation(S(:,:,1), [], [], [], 5, 'gaussian');
    figure; [PNCC2, lags] = plotSvPartialCrossCorrelation(S(:,:,2), [], [], [], 5, 'gaussian');
    PNCC1(:,:,6)=[]; PNCC2(:,:,6)=[];
    getCosSimilarity(PNCC1+U, PNCC2+U)
    
    load(['results/simsc/v09-' num2str(n) 'x2000-idx6-1-1-result.mat']);
    figure; [PNCC2, lags] = plotSvPartialCrossCorrelation(S(:,:,2), [], [], [], 5, 'gaussian');
    PNCC2(:,:,6)=[];
    getCosSimilarity(PNCC1+U, PNCC2+U)
    
    load(['results/simsc/mft-' num2str(n) 'x2000-idx6-1-1-result.mat']);
    figure; [PNCC2, lags] = plotSvPartialCrossCorrelation(S(:,:,2), [], [], [], 5, 'gaussian');
    PNCC2(:,:,6)=[];
    getCosSimilarity(PNCC1+U, PNCC2+U)
end

