%%
% get DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)

function wcI = calcDlcmWCI(netDLCM, nodeControl, inControl)
    wcI = calcDlcmWCIm1(netDLCM, nodeControl, inControl);
end
