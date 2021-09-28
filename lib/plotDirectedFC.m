%%
% plot directed FC matrix
% input:
%  dFC          directed FC matrix
%  name         directed FC name for graph title
%  range        plotting minimum and maximum range of dFC (default:0.5)

function plotDirectedFC(dFC, name, range)
    if nargin < 3
        range = 0.5;
    end
    % show directed FC
    if range <= 0
        sigma = std(dFC(:),1,'omitnan');
        avg = mean(dFC(:),'omitnan');
        dFC = (dFC - avg) / sigma;
        range = 3;
    end
    clims = [-range, range];
    imagesc(dFC,clims);
    daspect([1 1 1]);
    title([name ' directed FC Matrix']);
    xlabel('Source Nodes');
    ylabel('Target Nodes');
    colorbar;
end
