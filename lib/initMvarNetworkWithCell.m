%%
% calculate linear multivariate Vector Auto-Regression from Cells and Create mVAR network
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  uniqueDecimal   taking unique value from conjuncted time series (option)
%  verbose         show verbose log (default:false)

function net = initMvarNetworkWithCell(CX, CexSignal, nodeControl, exControl, lags, uniqueDecimal, verbose)
    if nargin < 7, verbose = false; end
    if nargin < 6, uniqueDecimal = 0; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, CexSignal = {}; end
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    sigLen = size(CX{1},2);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    inputNum = nodeNum + exNum;

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    b = cell(nodeNum,1);
    r = cell(nodeNum,1);
    T = cell(nodeNum,1);

    % set vector auto-regression (VAR) inputs
    idxs = cell(nodeNum,1);
    for n=1:nodeNum
        [~,idxs{n}] = find(control(n,:,:)==1);
    end
    allInLen = 0;
    for i=1:cxNum
        allInLen = allInLen + size(CX{i},2) - lags;
    end

    % calculate mean and covariance of each node
    Y = [];
    for i=1:cxNum, Y = [Y, CX{i}]; end
    cxM = mean(Y.');
    cxCov = cov(Y.');
    Y = []; % memory clear

    % apply the regress function
%    for n=1:nodeNum
    parfor n=1:nodeNum
        if verbose, disp(['calc node' num2str(n)]); end
        Xt = single(nan(allInLen,1));
        Xti = single(nan(allInLen,length(idxs{n})));
        xts = 1;

        % this is redundant for every node. but it is necessary to avoid
        % too much memory consumption
        for i=1:cxNum
            % set node input
            Y = single(CX{i});
            if exNum > 0
                Y = [Y; single(CexSignal{i})];
            end
            Y = flipud(Y.'); % need to flip signal

            sLen = size(Y,1);
            sl = sLen-lags;
            Yj = single(zeros(sl, lags*inputNum));
            for p=1:lags
                Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sl+p,:);
            end
            Xt(xts:xts+sl-1,:) = Y(1:sl,n);
            Xti(xts:xts+sl-1,:) = Yj(:,idxs{n});
            xts = xts + sl;
        end
        Y = [];  % clear memory
        Yj = []; % clear memory
        
        if uniqueDecimal > 0
            A = int32([Xt, Xti] / uniqueDecimal);
            A = unique(A,'rows');
            Xt = single(A(:,1)) * uniqueDecimal;
            Xti = single(A(:,2:end)) * uniqueDecimal;
            A = []; % clear memory
        end

        Xti = [Xti, ones(size(Xti,1),1)]; % add bias
        [b{n}, r{n}, T{n}] = regressLinear(Xt,Xti);
%{        
        [b2,~,r2,~] = regress(Xt,Xti);
        db2 = sum(abs(b{n}-b2));
        dr2 = sum(abs(r{n}-r2));
        disp(['node=' num2str(n) ', diffB=' num2str(db2) ', diffR=' num2str(dr2)]);
%}
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.sigLen = sigLen;
    net.cxM = cxM;
    net.cxCov = cxCov;
    net.lags = lags;
    net.bvec = b;
    net.rvec = r;
    net.Tvec = T;
end
