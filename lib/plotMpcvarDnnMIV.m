%%
% Plotting multivariate PC VAR DNN Mean Impact Value (MIV) matrix
% returns multivariate PC VAR DNN MIV matrix (MIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          trained multivariate PC VAR DNN network
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function MIV = plotMpcvarDnnMIV(X, exSignal, nodeControl, exControl, net, range, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, range = 10; end

    MIV = calcMpcvarDnnMIV(X, exSignal, nodeControl, exControl, net, isFullNode);
    if range <= 0
        sigma = std(MIV(:),1,'omitnan');
        avg = mean(MIV(:),'omitnan');
        MIV2 = (MIV - avg) / sigma;
        range = 3;
    else
        MIV2 = MIV;
    end
    clims = [-range, range];
    imagesc(MIV2,clims);
    daspect([1 1 1]);
    title('multivariate PC VAR DNN Mean Impact Value');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
