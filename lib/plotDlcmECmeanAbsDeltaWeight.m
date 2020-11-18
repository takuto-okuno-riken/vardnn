%%
% get DLCM effective connectivity matrix (dlEC) and its standard error matrix
% input:
%  netDLCM   trained DLCM network
%  range     plotting minimum and maximum range of GCI (default:0.25)

function [dlEC, dlECerr] = plotDlcmECmeanAbsDeltaWeight(netDLCM, range)
    if nargin < 2
        range = 0.2;
        rangemin = 0;
    end
    [dlEC, dlECerr] = getDlcmECmeanAbsDeltaWeight(netDLCM);
    % show effective conectivity of predicted node signals
    if range <= 0
        sigma = std(dlEC(:),1,'omitnan');
        avg = mean(dlEC(:),'omitnan');
        dlEC = (dlEC - avg) / sigma;
        range = 5;
        rangemin = -5;
    end
    clims = [rangemin,range];
    imagesc(dlEC,clims);
    daspect([1 1 1]);
    title('DL Effective Connectivity');
    colorbar;
end
