%%
% Plotting Time Shifted Correlation matrix
% returns Functional Connectivity (FC)
% input:
%  X      multivariate time series matrix (node x time series)
%  lags   number of lags for shifting (default:3)

function [tscFC] = plotTimeShiftedCorrelation_(X, lag)
    if nargin < 2
        lag = 3;
    end
    tscFC = calcTimeShiftedCorrelation_(X, lag);
    clims = [-1,1];
    imagesc(tscFC,clims);
    daspect([1 1 1]);
    title('Time Shifted Correlation (FC)');
    colorbar;
end
