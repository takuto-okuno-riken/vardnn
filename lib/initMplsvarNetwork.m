%%
% calculate linear multivariate PLS Vector Auto-Regression weights and Create mPLSVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)

function net = initMplsvarNetwork(X, exSignal, nodeControl, exControl, lags)
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

    b = cell(nodeNum,1);
    XL = cell(nodeNum,1);
    YL = cell(nodeNum,1);
    XS = cell(nodeNum,1);
    YS = cell(nodeNum,1);
    PCTVAR = cell(nodeNum,1);
    MSE = cell(nodeNum,1);
    stats = cell(nodeNum,1);

    % find component number
    ncomp = floor(inputNum * lags / 10);
    if ncomp < 5, ncomp = 5; end
    if ncomp > 50, ncomp = 50; end
    if ncomp > (sigLen-lags-1), ncomp = (sigLen-lags-1); end
    
    % first, calculate PLS vector auto-regression (VAR) without target
    Yj = zeros(sigLen-lags, lags*inputNum);
    for k=1:lags
        Yj(:,1+inputNum*(k-1):inputNum*k) = Y(1+k:sigLen-lags+k,:);
    end
    for i=1:nodeNum
        [~,idx] = find(control(i,:,:)==1);
        
        % PLS vector auto-regression (VAR)
        Xt = Y(1:sigLen-lags,i);
        Xti = Yj(:,idx);
        % apply the PLS regress function
        [XL{i},YL{i},XS{i},YS{i},b{i},PCTVAR{i},MSE{i},stats{i}] = plsregress(Xti,Xt,ncomp);
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.ncomp = ncomp;
    net.bvec = b;
    net.XL = XL;
    net.YL = YL;
    net.XS = XS;
    net.YS = YS;
    net.PCTVAR = PCTVAR;
    net.MSE = MSE;
    net.stats = stats;
end