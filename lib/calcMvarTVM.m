%%
% calc mVAR (multivaliate Vector Auto-Regression) T-value matrix
% returns 3D mVAR (multivaliate Vector Auto-Regression) T-value matrix (TV)(node x node x lags)
% input:
%  net          mVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [TV] = calcMvarTVM(net, nodeControl, exControl, isFullNode)
    if nargin < 4, isFullNode = 0; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    if ~isfield(net, 'Tvec')
        disp('Tvec is not found');
        return;
    end

    % get TV matrix
    nodeNum = net.nodeNum;
    inputNum = nodeNum + net.exNum;
    lags = net.lags;
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = inputNum; end

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, net.exNum, lags);

    TV = nan(nodeNum, inputNum, lags);
    for i=1:length(net.Tvec)
        T = net.Tvec{i};
        last = 1;
        for k=1:lags
            idx = find(control(i,:,k)==1);
            num = length(idx);
            TV(i,idx,k) = T(last:last+num-1);
            last = last+num;
        end
    end
    TV = TV(:,1:nodeMax,:);
end
