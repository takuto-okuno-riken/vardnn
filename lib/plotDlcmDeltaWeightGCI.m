%%
% plot DLCM delta weight granger causality index matrix
% input:
%  netDLCM   trained DLCM network
%  range     plotting minimum and maximum range of GCI (default:10)

function [gcI] = plotDlcmDeltaWeightGCI(netDLCM, range)
    if nargin < 2
        range = 0.1;
    end
    gcI = getDlcmDeltaWeightGCI(netDLCM);
    % show DLCM granger causality of predicted node signals
    if range <= 0
        amax = abs(max(max(gcI)));
        amin = abs(min(min(gcI)));
        range = max(amax,amin);
    end
    clims = [-range,range];
    imagesc(gcI,clims);
    daspect([1 1 1]);
    title('DLCM delta weight Granger Causality');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
