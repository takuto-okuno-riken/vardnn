%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM   trained DLCM network
%  range     plotting minimum and maximum range of GCI (default:0.25)

function [dlEC, dlECerr] = plotDlcmECmeanWeight(netDLCM, range)
    if nargin < 2
        range = 0.25;
    end
    [dlEC, dlECerr] = getDlcmECmeanWeight(netDLCM);
    % show effective conectivity of predicted node signals
    if range <= 0
        amax = abs(max(max(dlEC)));
        amin = abs(min(min(dlEC)));
        range = max(amax,amin);
    end
    clims = [-range,range];
    imagesc(dlEC,clims);
    daspect([1 1 1]);
    title('DL Effective Connectivity');
    colorbar;
end
