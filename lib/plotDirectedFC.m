%%
% plot directed FC matrix
% input:
%  dFC          directed FC matrix
%  name         EC name for graph title
%  range        plotting minimum and maximum range of dFC (default:0.5)
%  rowcut       cut bottom rows of result CI matrix (default:0)

function plotDirectedFC(dFC, name, range, rowcut)
    if nargin < 4
        rowcut = 0;
    end
    if nargin < 3
        range = 0.5;
    end
    % show effective connectivity
    if range <= 0
        sigma = std(dFC(:),1,'omitnan');
        avg = mean(dFC(:),'omitnan');
        dFC = (dFC - avg) / sigma;
        range = 3;
    end
    if rowcut>0, dFC(end-rowcut+1:end,:) = []; end
    clims = [-range, range];
    imagesc(dFC,clims);
    daspect([1 1 1]);
    title([name ' directed FC Matrix']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
