%%
% plot DLCM weight granger causality index matrix
% input:
%  netDLCM   trained DLCM network
%  range     plotting minimum and maximum range of GCI (default:10)

function [gcI] = plotDlcmWeightGCI(netDLCM, range)
    if nargin < 2
        range = 0.1;
    end
    gcI = getDlcmWeightGCI(netDLCM);
    % show DLCM granger causality of predicted node signals
    if range <= 0
        amax = abs(max(max(gcI)));
        amin = abs(min(min(gcI)));
        range = max(amax,amin);
    end
    clims = [-range,range];
    imagesc(gcI,clims);
    daspect([1 1 1]);
    title('DLCM weight Granger Causality');
    colorbar;
end
