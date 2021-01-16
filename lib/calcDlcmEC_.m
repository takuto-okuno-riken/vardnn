%%
% calc DLCM effective connectivity matrix (EC) and impaired node signals (ecNS)
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [EC, ecNS] = calcDlcmEC_(netDLCM, nodeControl, inControl, isFullNode)
    [EC, ecNS] = calcDlcmWCIdm123a(netDLCM, nodeControl, inControl, isFullNode);
end
