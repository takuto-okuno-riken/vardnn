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
    p = lags;

    % set node input
    Y = [X; exSignal];
    nodeMax = nodeNum + exNum;
    
    b = cell(nodeNum,nodeMax);
    XL = cell(nodeNum,nodeMax);
    YL = cell(nodeNum,nodeMax);
    XS = cell(nodeNum,nodeMax);
    YS = cell(nodeNum,nodeMax);
    PCTVAR = cell(nodeNum,nodeMax);
    MSE = cell(nodeNum,nodeMax);
    stats = cell(nodeNum,nodeMax);

    % find component number
    ncomp = floor(2 * p / 2);
    if ncomp < 2, ncomp = 2; end
    if ncomp > 50, ncomp = 50; end

    for i=1:nodeNum
        for j=1:nodeMax
            if i==j, continue; end
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            Y1 = flipud(Y(i,:));
            Y2 = flipud(Y(j,:));

            % autoregression plus other regression
            Yt = Y2(1:sigLen-p).'; % TODO: X1 & X2 opposite ??
            Yti = ones(sigLen-p, p*2);
            for k=1:p
                Yti(:,k) = Y2(k+1:sigLen-p+k);
                Yti(:,p+k) = Y1(k+1:sigLen-p+k);
            end
            % apply the PLS regress function
            [XL{i,j},YL{i,j},XS{i,j},YS{i,j},b{i,j},PCTVAR{i,j},MSE{i,j},stats{i,j}] = plsregress(Yti,Yt,ncomp);
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