%%
% plot DLCM weight causality index matrix
% input:
%  netDLCM      trained DLCM network
%  nodeControl  node control matrix (node x node) (optional)
%  inControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  range        plotting minimum and maximum range of GCI (default:0.5)
%  rowcut       cut bottom rows of result CI matrix (default:0)

function [wcI] = plotDlcmWCI(netDLCM, nodeControl, inControl, range, rowcut)
    if nargin < 5
        rowcut = 0;
    end
    if nargin < 4
        range = 0.5;
    end
    if nargin < 3
        inControl = [];
    end
    if nargin < 2
        nodeControl = [];
    end
    nodeNum = length(netDLCM.nodeNetwork);
    wcI = calcDlcmWCI(netDLCM, nodeControl, inControl);
    % show DLCM weight causality of predicted node signals
    if range <= 0
        sigma = std(wcI(:),'omitnan');
        avg = mean(wcI(:),'omitnan');
        wcI = (wcI - avg) / sigma;
        range = 3;
    end
    wcI = wcI(:, 1:nodeNum); % TODO: probably, calcDlcmWeightCI should be changed.
    if rowcut>0, wcI(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(wcI,clims);
    daspect([1 1 1]);
    title('DLCM Weight Causality Index');
    colorbar;
end
