%%
% plot Effective Connectivity matrix
% input:
%  EC           Effective Connectivity matrix
%  name         EC name for graph title
%  range        plotting minimum and maximum range of GCI (default:0.5)
%  rowcut       cut bottom rows of result CI matrix (default:0)

function plotEC(EC, name, range, rowcut)
    if nargin < 4
        rowcut = 0;
    end
    if nargin < 3
        range = 0.5;
    end
    % show effective connectivity
    if range <= 0
        sigma = std(EC(:),1,'omitnan');
        avg = mean(EC(:),'omitnan');
        EC = (EC - avg) / sigma;
        range = 3;
    end
    if rowcut>0, EC(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(EC,clims);
    daspect([1 1 1]);
    title([name ' Effective Connectivity Matrix']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
