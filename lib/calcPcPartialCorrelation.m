%%
% Caluclate Principal Component Partial Correlation
% returns Principal Component Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  explainedTh  explained threshold for PCA components (default:0.99)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [PC] = calcPcPartialCorrelation(X, exSignal, nodeControl, exControl, explainedTh, isFullNode)
    if nargin < 6, isFullNode = 0; end
    if nargin < 5, explainedTh = 0.99; end
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
    for i=1:nodeMax
        A = unique(Y(i,:));
        if length(A)==1
            Y(i,mod(i,sigLen)+1) = A + 1.0e-8;
        end
    end

    fullIdx = 1:nodeMax;
    PC = nan(nodeNum,nodeMax);
    parfor i=1:nodeNum
        if ~isempty(nodeControl), nidx = find(nodeControl(i,:)==0); else nidx = []; end
        if ~isempty(exControl), eidx = find(exControl(i,:)==0); else eidx = []; end
        if ~isempty(eidx), eidx = eidx + nodeNum; end
        nodeIdx = setdiff(fullIdx,[nidx, eidx, i]);
        x = Y(i,:).';
        
        A = nan(1,nodeMax,class(X));
        for j=i:nodeMax
            if j<=nodeNum && ~isempty(nodeControl) && nodeControl(i,j) == 0, continue; end
            if j>nodeNum && ~isempty(exControl) && exControl(i,j-nodeNum) == 0, continue; end
            
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
            [~, ~, perm, RiQ] = regressPrepare(pcXti);
            [~, r1] = regressLinear(x, pcXti, [], [], perm, RiQ);
            [~, r2] = regressLinear(Y(j,:).', pcXti, [], [], perm, RiQ);
            
            A(j) = (r1.'*r2) / (sqrt(r1.'*r1)*sqrt(r2.'*r2));
        end
        PC(i,:) = A;
    end
    A = ones(nodeNum,'logical'); A = tril(A,-1);
    idx = find(A==1);
    B = PC(:,1:nodeNum); C = B';
    B(idx) = C(idx);
    PC(:,1:nodeNum) = B;
    
    % output control
    PC = PC(1:nodeNum,:);
    if isFullNode == 0
        PC = PC(:,1:nodeNum);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        PC(:,1:nodeNum) = PC(:,1:nodeNum) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        PC(:,nodeNum+1:end) = PC(:,nodeNum+1:end) .* exControl;
    end
end
