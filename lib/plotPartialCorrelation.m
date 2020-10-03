%%
% Plot Partial Correlation
% returns Partial Correlation (PC) and p-values (P)
% input:
%  X       multivariate time series matrix (node x time series)
%  rowcut  cut bottom rows of result gcI matris (default:0)

function [PC, P] = plotPartialCorrelation(X, rowcut)
    if nargin < 2
        rowcut = 0;
    end
    % show partial correlation
    [PC, P] = calcPartialCorrelation(X);
    if rowcut>0, PC(end-rowcut+1:end,:) = []; end
    clims = [-1,1];
    imagesc(PC,clims);
    daspect([1 1 1]);
    title('Partial Correlation');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
