%%
% Plotting mVAR (multivaliate Vector Auto-Regression) T-value matrix
% returns 3D mVAR (multivaliate Vector Auto-Regression) T-value matrix (TV)(node x node x lags)
% input:
%  net          mVAR network
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  range        plotting minimum and maximum range of T-value (default:20)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (optional)

function [TV] = plotMvarTVM(net, nodeControl, exControl, range, isFullNode)
    if nargin < 5, isFullNode = 0; end
    if nargin < 4, range = 20; end
    if nargin < 3, exControl = []; end
    if nargin < 2, nodeControl = []; end

    TV = calcMvarTVM(net, nodeControl, exControl, isFullNode);

    if range <= 0
        sigma = std(TV(:),1,'omitnan');
        avg = mean(TV(:),'omitnan');
        DI2 = (TV - avg) / sigma;
        range = 3;
    else
        DI2 = TV;
    end
    n = ceil(sqrt(net.lags));
    for i=1:net.lags
        subplot(n,n,i)
        clims = [-range, range];
        imagesc(DI2(:,:,i),clims);
        daspect([1 1 1]);
        title(['mVAR T-value matrix (' num2str(i) ')']);
        xlabel('Source Nodes');
        ylabel('Target Nodes');
    end
    colorbar;
end
