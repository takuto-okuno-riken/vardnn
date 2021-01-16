%%
% calc DLCM effective connectivity matrix (EC) and impaired node signals (ecNS)
% input:
%  netDLCM      trained DLCM network

function [EC, ecNS] = calcDlcmEC(netDLCM)
    [EC, ecNS] = calcDlcmEC_(netDLCM, [], []);
end
