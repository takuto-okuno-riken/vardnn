%%
% Caluclate Partial Correlation by Inverse Covariance
% returns Partial Correlation (PC)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)
%  usegpu       use gpu calculation (default:false)

function PC = calcPartialCorrelation_(X, exSignal, nodeControl, exControl, isFullNode, usegpu)
    if nargin < 6, usegpu = false; end
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end
    nodeNum = size(X,1);
    exNum = size(exSignal,1);
    nodeMax = nodeNum + exNum;

    % set node input
    Y = [X; exSignal];
    if isa(Y,'half')
        Y = double(Y); % half cannot be calculated by qr()
    end
    if usegpu
        Y = gpuArray(single(Y));
    end

    % using matrix inversion
    PC = nan(nodeMax,nodeMax,class(X));
    P = zeros(nodeMax,nodeMax);
    C = cov(Y',1);
    [sz1,sz2] = size(C);
    [Q,R,perm] = qr(C,0);
    p = sum(abs(diag(R)) > max(sz1,sz2)*eps(R(1))); % 2 steps differ than zero
    if p < sz2
        R = R(1:p,1:p);
        Q = Q(:,1:p);
        perm = perm(1:p);
    end
    Ci = inv(R) * Q'; % get precision matrix
    P(perm,:) = Ci;
    clear Q; clear R; clear perm; 
    Dp = diag(P);
    pii = repmat(Dp(:),1,nodeMax);
    pjj = repmat(Dp(:)',nodeMax,1);
    P2 = -P ./ sqrt(pii .* pjj);
    clear pii; clear pjj;
    P2(eye(nodeMax)==1) = 1; % replace diag elements
    PC(:,:) = P2;

    % output control
    PC = PC(1:nodeNum,:,:);
    if isFullNode == 0
        PC = PC(:,1:nodeNum,:);
    end
    if ~isempty(nodeControl)
        nodeControl=double(nodeControl); nodeControl(nodeControl==0) = nan;
        PC(:,1:nodeNum,:) = PC(:,1:nodeNum,:) .* nodeControl;
    end
    if ~isempty(exControl) && isFullNode > 0
        exControl=double(exControl); exControl(exControl==0) = nan;
        PC(:,nodeNum+1:end,:) = PC(:,nodeNum+1:end,:) .* exControl;
    end
end
