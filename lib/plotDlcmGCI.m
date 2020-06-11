%%
% Plotting DLCM Granger causality Index matrix
% returns DLCM Granger Causality Index (gcI)
% input:
%  X          multivariate time series matrix (node x time series)
%  inSignal   multivariate time series matrix (exogenous input x time series) (optional)
%  inControl  exogenous input control matrix for each node (node x exogenous input) (optional)
%  range   plotting minimum and maximum range of GCI (default:10)
%  rowcut  cut bottom rows of result gcI matris (default:0)

function [gcI] = plotDlcmGCI(X, inSignal, inControl, netDLCM, range, rowcut)
    if nargin < 6
        rowcut = 0;
    end
    if nargin < 5
        range = 10;
    end
    gcI = calcDlcmGCI(X, inSignal, inControl, netDLCM);
    if rowcut>0, gcI(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(gcI,clims);
    daspect([1 1 1]);
    title('DLCM Granger Causality');
    colorbar;
end
