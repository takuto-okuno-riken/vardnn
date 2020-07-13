%%
% Plot DCM Effective Connectivity
% returns DCM Effective Connectivity (EC)
% input:
%  A       DCM connection weight matrix A
%  range   plotting minimum and maximum range of GCI (default:1)

function plotDcmEC(A, range)
    if nargin < 2
        range = 1;
    end 
    if range <= 0
        amax = abs(max(max(A)));
        amin = abs(min(min(A)));
        range = max(amax,amin);
    end
    clims = [-range, range];
    imagesc(A,clims);
    daspect([1 1 1]);
    title('DCM Effective Connectivity');
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
