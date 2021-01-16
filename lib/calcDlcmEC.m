%%
% calc DLCM effective connectivity matrix (EC) and impaired node signals (ecNS)
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [EC, ecNS] = calcDlcmEC(netDLCM, nodeControl, exControl, isFullNode)
    if nargin < 4
        isFullNode = 0;
    end
    if nargin < 3
        exControl = [];
    end
    if nargin < 2
        nodeControl = [];
    end
    [EC, ecNS] = calcDlcmWCIdm123a(netDLCM, nodeControl, exControl, isFullNode);
end
