%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM   trained DLCM network
%  range     plotting minimum and maximum range of GCI (default:0.25)

function [dlEC, dlECerr] = plotDlcmECmeanAbsDeltaWeight(netDLCM, range)
    if nargin < 2
        range = 0.2;
    end
    [dlEC, dlECerr] = getDlcmECmeanAbsDeltaWeight(netDLCM);
    % show effective conectivity of predicted node signals
    if range <= 0
        amax = abs(max(max(dlEC)));
        amin = abs(min(min(dlEC)));
        range = max(amax,amin);
    end
    clims = [0,range];
    imagesc(dlEC,clims);
    daspect([1 1 1]);
    title('DL Effective Connectivity');
    colorbar;
end
