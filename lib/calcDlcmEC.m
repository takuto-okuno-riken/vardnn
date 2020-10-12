%%
% calc DLCM effective connectivity matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function EC = calcDlcmEC(netDLCM, nodeControl, inControl)
    EC = calcDlcmWCIdm123a(netDLCM, nodeControl, inControl);
end
