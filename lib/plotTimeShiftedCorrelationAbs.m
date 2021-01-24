%%
% Plotting Time Shifted Correlation (Abs) matrix
% returns Functional Connectivity (FC) and p-values (P)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  lags         number of lags for time shift (default:3)
%  range        plotting minimum and maximum range of FC (default:1)
%               if range==0, range shows standard deviation [-3 sigma, 3 sigma]
%  isFullNode   return both node & exogenous causality matrix (default:0)

function [tscFC, P] = plotTimeShiftedCorrelationAbs(X, exSignal, nodeControl, exControl, lags, range, isFullNode)
    if nargin < 7, isFullNode = 0; end
    if nargin < 6, range = 1; end
    if nargin < 5, lags = 3; end
    if nargin < 4, exControl = []; end
    if nargin < 3, nodeControl = []; end
    if nargin < 2, exSignal = []; end

    [tscFC, P] = calcTimeShiftedCorrelationAbs(X, exSignal, nodeControl, exControl, lags, isFullNode);

    if range <= 0
        sigma = std(tscFC(:),1,'omitnan');
        avg = mean(tscFC(:),'omitnan');
        tscFC2 = (tscFC - avg) / sigma;
        range = 3;
    else
        tscFC2 = tscFC;
    end
    clims = [-range, range];
    imagesc(tscFC2,clims);
    daspect([1 1 1]);
    title('Time Shifted Correlation (Abs) (FC)');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
