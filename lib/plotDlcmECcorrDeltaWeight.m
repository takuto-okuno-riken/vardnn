%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM   trained DLCM network
%  range     plotting minimum and maximum range of GCI (default:0.25)

function [dlEC] = plotDlcmECcorrDeltaWeight(netDLCM)
    [dlEC] = getDlcmECcorrDeltaWeight(netDLCM);
    % show effective conectivity of predicted node signals
    clims = [-1,1];
    imagesc(dlEC,clims);
    daspect([1 1 1]);
    title('DL Effective Connectivity');
    colorbar;
end
