
function testMrfvarDI3
    % load signals
    load('test/testTrain-rand500-uniform.mat');
    siOrg = si;
    nodeNum = 8;
    exNum = 0;
    sigLen = 200;
    lags = 1;
    si = siOrg(1:nodeNum,1:sigLen);
    exSignal = [];
    % control is all positive input
    nodeControl = [];
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1); % lag=1
    si(4,3:end) = si(6,1:sigLen-2); % lag=2

    %% test pattern 1 
    % init multivariate RF VAR network
    net = initMrfvarNetwork(si, exSignal, nodeControl, exControl, lags);

    % show mRFVAR-MIV3
    [MIV3,MAIV3] = calcMrfvarMIV3(si, exSignal, nodeControl, exControl, net, 0);
    figure; plotDirectedFC3(MIV3,'mRFVAR-MIV',0);
    figure; plotDirectedFC3(MAIV3,'mRFVAR-MAIV',0);

    %% test pattern 2 -- exogenous signals
    % do training or load multivariate RF VAR network
    lags = 3;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum,lags);
    si(3,2:end) = exSignal(1,1:sigLen-1); % lag=1
    si(5,4:end) = exSignal(2,1:sigLen-3); % lag=3
    % init multivariate RF VAR network
    net = initMrfvarNetwork(si, exSignal, nodeControl, exControl, lags);

    % show mRFVAR-MIV3
    [MIV3,MAIV3] = calcMrfvarMIV3(si, exSignal, nodeControl, exControl, net, 1);
    figure; plotDirectedFC3(MIV3,'mRFVAR-MIV',0);
    figure; plotDirectedFC3(MAIV3,'mRFVAR-MAIV',0);

    %% test pattern 3 -- input value selection test
    % do training or load multivariate RF VAR network
    lags = 10;
    exNum = 2;
    exSignal = siOrg(nodeNum+1:nodeNum+exNum,1:sigLen);
    % control is all positive input
    exControl = ones(nodeNum,exNum,lags);
    si(3,4:end) = exSignal(1,1:sigLen-3); % lag=3
    si(5,8:end) = exSignal(2,1:sigLen-7); % lag=7

    % init multivariate RF VAR network
    net = initMrfvarNetwork(si, exSignal, nodeControl, exControl, lags);

    % input value selection by RFVAR-MIV
    ivsName = 'ivs2';
    rate = 8 / (10*10*lags);
    [MIV3,MAIV3] = calcMrfvarMIV3(si, exSignal, [], [], net, 1);
    control3 = inputValueSelectionFromDFC3(MIV3, rate);
    nodeControl = control3(:,1:nodeNum,:);
    exControl = control3(:,nodeNum+1:nodeNum+exNum,:);

    % init multivariate RF VAR network
    net2 = initMrfvarNetwork(si, exSignal, nodeControl, exControl, lags);
    [MIV3,MAIV3] = calcMrfvarMIV3(si, exSignal, nodeControl, exControl, net2, 1);
    figure; plotDirectedFC3(MIV3,[ivsName ' mRFVAR-MIV'],0);
    figure; plotDirectedFC3(MAIV3,[ivsName ' mRFVAR-MAIV'],0);

    % input value selection by RFVAR-DI
    ivsName = 'ivs2';
    [MIV3,MAIV3] = calcMrfvarMIV3(si, exSignal, [], [], net, 1);
    control3 = inputValueSelectionFromDFC3(MAIV3, rate);
    nodeControl = control3(:,1:nodeNum,:);
    exControl = control3(:,nodeNum+1:nodeNum+exNum,:);

    % init multivariate RF VAR network
    net2 = initMrfvarNetwork(si, exSignal, nodeControl, exControl, lags);
    [MIV3,MAIV3] = calcMrfvarMIV3(si, exSignal, nodeControl, exControl, net2, 1);
    figure; plotDirectedFC3(MIV3,[ivsName ' mRFVAR-MIV'],0);
    figure; plotDirectedFC3(MAIV3,[ivsName ' mRFVAR-MAIV'],0);
end



