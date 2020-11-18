%%
% calc DLCM effective connectivity matrix (EC) and impaired node signals (ecNS)
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function [EC, ecNS] = calcDlcmEC(netDLCM, nodeControl, inControl)
    [EC, ecNS] = calcDlcmWCIdm123a(netDLCM, nodeControl, inControl);
end
