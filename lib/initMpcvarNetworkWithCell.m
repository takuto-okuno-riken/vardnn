%%
% calculate linear multivariate Principal Component Vector Auto-Regression weights and Create mPCVAR network
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)
%  explainedTh     explained threshold for PCA components (default:0.99)

function net = initMpcvarNetworkWithCell(CX, CexSignal, nodeControl, exControl, lags, explainedTh)
    if nargin < 6, explainedTh = 0.99; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, CexSignal = {}; end
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    inputNum = nodeNum + exNum;
    expTh = explainedTh * 100;

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    coeff = cell(nodeNum,1);
%    score = cell(nodeNum,1);
    latent = cell(nodeNum,1);
    explained = cell(nodeNum,1);
    mu = cell(nodeNum,1);
    maxComp = cell(nodeNum,1);
    b = cell(nodeNum,1);
    r = cell(nodeNum,1);
    stats = cell(nodeNum,1);
    
    % init response variables
    for n=1:nodeNum, Xt{n} = []; Xti{n} = []; end

    % set vector auto-regression (VAR) inputs
    for i=1:cxNum
        % set node input
        Y = single(CX{i});
        if exNum > 0
            Y = [Y; single(CexSignal{i})];
        end
        Y = flipud(Y.'); % need to flip signal

        sLen = size(Y,1);
        Yj = single(zeros(sLen-lags, lags*inputNum));
        for p=1:lags
            Yj(:,1+inputNum*(p-1):inputNum*p) = Y(1+p:sLen-lags+p,:);
        end
        for n=1:nodeNum
            [~,idx] = find(control(n,:,:)==1);
            Xt{n} = [Xt{n}; Y(1:sLen-lags,n)];
            Xti{n} = cat(1,Xti{n},Yj(:,idx));
        end
    end

    % apply the Principal Component Regress function
    for n=1:nodeNum
        [coeff{n},score,latent{n},~,explained{n},mu{n}] = pca(Xti{n}); % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);

        % find 99% component range
        expTotal = 0;
        maxComp{n} = size(score,2);
        for j=1:size(Xti{n},2)
            expTotal = expTotal + explained{n}(j);
            if expTotal >= expTh
                maxComp{n} = j;
                break;
            end
        end
        pcXti = [score(:,1:maxComp{n}), ones(size(Xti{n},1),1)]; % might not be good to add bias
        [b{n},~,r{n},~,stats{n}] = regress(Xt{n}, pcXti);
        b{n} = single(b{n});
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.coeff = coeff;
    net.latent = latent;
    net.explained = explained;
    net.explainedTh = explainedTh;
    net.mu = mu;
    net.maxComp = maxComp;
    net.bvec = b;
    net.rvec = r;
    net.stats = stats;
end