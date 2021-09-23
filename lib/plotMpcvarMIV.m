%%
% Plotting mPCVAR (multivaliate Principal Component Vector Auto-Regression) Mean Inpact Value (MIV)
% returns Mean Inpact Value matrix (MIV)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mPCVAR network
%  range        plotting minimum and maximum range of GCI (default:10)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function MIV = plotMpcvarMIV(X, exSignal, nodeControl, exControl, net, range, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, range = 10; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    MIV = calcMpcvarMIV(X, exSignal, nodeControl, exControl, net, isFullNode);
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
    title('mPCVAR Mean Inpact Value');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
