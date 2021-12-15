%%
% calculate linear multivariate Vector Auto-Regression from Cells and Create mVAR network
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  CexSignal       cells of multivariate time series matrix {exogenous input x time series} (default:{})
%  nodeControl     node control matrix (node x node) (default:[])
%  exControl       exogenous input control matrix for each node (node x exogenous input) (default:[])
%  lags            number of lags for autoregression (default:3)

function net = initMvarNetworkWithCell(CX, CexSignal, nodeControl, exControl, lags)
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

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

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
    
    % apply the regress function
    for n=1:nodeNum
        Xti{n} = [Xti{n}, ones(size(Xti{n},1),1)]; % add bias
        [b{n},~,r{n},~,stats{n}] = regress(Xt{n},Xti{n});
        b{n} = single(b{n});
    end
    net.nodeNum = nodeNum;
    net.exNum = exNum;
    net.lags = lags;
    net.bvec = b;
    net.rvec = r;
    net.stats = stats;
end