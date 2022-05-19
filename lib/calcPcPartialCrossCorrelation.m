%%
% Caluclate Principal Component normalized Partial Cross-Correlation (PcPCC)
% returns Principal Component normalized Partial Cross-Correlation (PCC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  maxlag       maxlag of normalized cross-correlation [-maxlag, maxlag] (default:5)
%  explainedTh  explained threshold for PCA components (default:0.99)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [NPCC, lags] = calcPcPartialCrossCorrelation(X, exSignal, nodeControl, exControl, maxlag, explainedTh, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, explainedTh = 0.99; end
    if nargin < 5, maxlag = 5; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;
    expTh = explainedTh * 100;

    % set node input
    Y = [X; exSignal];

    % check all same value or not
    uLen = zeros(nodeMax,1);
    for i=1:nodeMax
        uLen(i) = length(unique(Y(i,:)));
    end
    
    NPCC = nan(nodeNum,nodeMax,maxlag*2+1);
    fullIdx = 1:nodeMax;
%    for i=1:nodeNum
    parfor i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);

        A = nan(nodeMax,maxlag*2+1);
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end

            x = Y(i,:).';
            y = Y(j,:).';
            idx = setdiff(nodeIdx,j);
            z = Y(idx,:).';

            % apply the Principal Component Regress function
            [coeff,score,latent,~,explained,mu] = pca(z); % relation : Xti == score * coeff.' + repmat(mu,size(score,1),1);
            % find 99% component range
            expTotal = 0;
            maxComp = size(score,2);
            for k=1:size(z,2)
                expTotal = expTotal + explained(k);
                if expTotal >= expTh
                    maxComp = k;
                    break;
                end
            end
            pcXti = [score(:,1:maxComp), ones(sigLen,1)]; % might not be good to add bias
            [Q, R, perm, RiQ] = regressPrepare(pcXti);
            [~, r1] = regressLinear(x, pcXti, Q, R, perm, RiQ);
            [~, r2] = regressLinear(y, pcXti, Q, R, perm, RiQ);

            [A(j,:), ~] = xcov(r1,r2,maxlag,'normalized');
        end
        NPCC(i,:,:) = A;
    end
    for i=1:nodeNum
        for j=i:nodeMax, NPCC(j,i,:) = NPCC(i,j,:); end
    end
    lags = -maxlag:maxlag;

    % output control
    NPCC = NPCC(1:nodeNum,:,:);
    if isFullNode == 0
        NPCC = NPCC(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        NPCC(:,1:nodeNum,:) = NPCC(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        NPCC(:,nodeNum+1:end,:) = NPCC(:,nodeNum+1:end,:) .* exControl;
    end
end
