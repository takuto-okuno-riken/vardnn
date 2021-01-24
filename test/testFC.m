function testFC
    nodeNum = 8;  % node number
    exNum = 2;    % exogenous number
    sigLen = 100; % signal length
    lags = 2;

    % generate random signals
    X = rand(nodeNum, sigLen); 
    exSignal = rand(exNum, sigLen);

    % copy signal 6->2, 6->4, 7 invert 3, ex1->1, ex2==5,
    X(2,2:end) = X(6,1:end-1);
    X(4,2:end) = X(6,1:end-1);
    X(3,:) = 1 - X(7,:);
    X(1,2:end) = exSignal(1,1:end-1);
    X(5,:) = exSignal(2,:);

    % plot (no exogenous)
    figure; FC = plotFunctionalConnectivity(X); % calc FC
    figure; FCa = plotFunctionalConnectivityAbs(X); % calc FC (Abs)
    figure; tsCr = plotTimeShiftedCorrelation(X, [], [], [], lags);
    figure; tsCra = plotTimeShiftedCorrelationAbs(X, [], [], [], lags);

    % plot with exogenous
    figure; FC = plotFunctionalConnectivity(X, exSignal, [], [], 1); % calc FC
    figure; FCa = plotFunctionalConnectivityAbs(X, exSignal, [], [], 1); % calc FC (Abs)
    figure; tsCr = plotTimeShiftedCorrelation(X, exSignal, [], [], lags, 0, 1);
    figure; tsCra = plotTimeShiftedCorrelationAbs(X, exSignal, [], [], lags, 0, 1);
end