%%
% plot 3D directed FC matrix
% input:
%  dFC          3D directed FC matrix
%  name         directed FC name for graph title
%  range        plotting minimum and maximum range of dFC (default:0.5)

function plotDirectedFC3(dFC, name, range)
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

    lags = size(dFC,3);
    n = ceil(sqrt(lags));
    for i=1:lags
        subplot(n,n,i)
        clims = [-range, range];
        imagesc(dFC(:,:,i),clims);
        daspect([1 1 1]);
        title([name '(' num2str(i) ')']);
        xlabel('Source Nodes');
        ylabel('Target Nodes');
    end
    colorbar;
end
