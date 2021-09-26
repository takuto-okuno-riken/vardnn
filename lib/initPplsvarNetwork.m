%%
% calculate linear pairwise PLS Vector Auto-Regression weights and Create pPLSVAR network
% input:
%  X               multivariate time series matrix (node x time series)
%  exSignal        multivariate time series matrix (exogenous input x time series) (default:[])
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)

function net = initPplsvarNetwork(X, exSignal, nodeControl, exControl, lags, explainedTh)
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    Y = flipud(Y.'); % need to flip signal

    % set control 3D matrix (node x node x lags)
    [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    b = cell(nodeNum,nodeMax);
    XL = cell(nodeNum,nodeMax);
    YL = cell(nodeNum,nodeMax);
    XS = cell(nodeNum,nodeMax);
    YS = cell(nodeNum,nodeMax);
    PCTVAR = cell(nodeNum,nodeMax);
    MSE = cell(nodeNum,nodeMax);
    stats = cell(nodeNum,nodeMax);

    for i=1:nodeNum
        [~,idx] = find(control(i,i,:)==1);
        Xt = Y(1:sigLen-lags,i);
        Yi = zeros(sigLen-lags, lags);
        for k=1:lags, Yi(:,k) = Y(1+k:sigLen-lags+k,i); end
        Xti = Yi(:,idx);

        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~any(nodeControl(i,j,:),'all'), continue; end
            if j>nodeNum && ~any(exControl(i,j-nodeNum,:),'all'), continue; end

            % autoregression plus other regression
            [~,idx] = find(control(i,j,:)==1);
            Yj = zeros(sigLen-lags, lags);
            for k=1:lags, Yj(:,k) = Y(1+k:sigLen-lags+k,j); end
            Xtj = [Xti, Yj(:,idx)];

            % find component number
            ncomp = floor(size(Xtj,2) / 2);
            if ncomp < 2, ncomp = 2; end
            if ncomp > 50, ncomp = 50; end

            % apply the PLS regress function
            [XL{i,j},YL{i,j},XS{i,j},YS{i,j},b{i,j},PCTVAR{i,j},MSE{i,j},stats{i,j}] = plsregress(Xtj,Xt,ncomp);
        end
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