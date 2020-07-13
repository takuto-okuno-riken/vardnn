%%
% Plot Functional Connectivity
% returns Functional Connectivity (FC)
% input:
%  X       multivariate time series matrix (node x time series)
%  rowcut  cut bottom rows of result gcI matris (default:0)

function [FC] = plotFunctionalConnectivity(X, rowcut)
    if nargin < 2
        rowcut = 0;
    end
    % show functional conectivity
    FC = calcFunctionalConnectivity(X);
    if rowcut>0, FC(end-rowcut+1:end,:) = []; end
    clims = [-1,1];
    imagesc(FC,clims);
    daspect([1 1 1]);
    title('Functional Connectivity');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
