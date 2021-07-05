% Before using this function, download Granger-causality codes from
% https://github.com/USC-Melady/Granger-causality
% and add a path "Granger-causality-master" and sub folders. And also download glmnet_matlab code from
% https://github.com/growlix/glmnet_matlab
% and add a path "glmnet_matlab-master" and sub folders. (Original glmnet was for windows7 and mex do not work)

function testLassoGranger
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
    exControl = [];
    si(2,2:end) = si(6,1:sigLen-1);
    si(4,2:end) = si(6,1:sigLen-1); % caution! node 2 & 4 is Multicollinearity case (correlated)

    %% test pattern 1 
    lambda = 1e-2;
    N = nodeNum;
    cause = zeros(N, N);
    for in = 1:N
        index = [in, 1:(in-1), (in+1):N];
        [~, temp] = lassoGranger(si(index, :), lags, lambda, 'l');
        cause(in, :) = temp([2:in, 1, (in+1):N])';
    end
    figure; plotDirectedFC(cause, 'Lasso Granger', 0);
end

