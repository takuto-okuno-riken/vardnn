%%
% calculate multivariate SVM(support vector machine) Vector Auto-Regression weights and Create mSvmVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  kernel          kernel for SVM (default:'linear', 'gaussian', 'rbf')
%  kernelScale     kernelScale for SVM (default:'auto', 1)

function net = initMsvmvarNetwork(X, exSignal, nodeControl, exControl, lags, kernel, kernelScale)
    if nargin < 7, kernelScale = 'auto'; end
    if nargin < 6, kernel = 'linear'; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    mdl = cell(nodeNum,1);

    % first, calculate SVM vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);

        % vector auto-regression (SVM VAR)
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);
        % apply the regress function
        mdl{i} = fitrsvm(Xti,Xt,'KernelFunction',kernel,'KernelScale',kernelScale); %,'Standardize',true); % bias will be calcurated
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.mdl = mdl;
    net.kernel = kernel;
    net.kernelScale = kernelScale;
end