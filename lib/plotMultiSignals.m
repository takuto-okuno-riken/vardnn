%%
% Prot multi signals
% input:
%  X          3D multivariate time series matrix (node x time series x group)
%  showEach   showSubplot (0) or show each graph (1) (default:0)
%  yRange     ylim range (default:[-0.2 1])

function plotMultiSignals(X, showEach, yRange)
    if nargin < 4, yRange = [-0.2 1]; end
    if nargin < 3, showEach = 0; end

    nodeNum = size(X,1);
    gNum = size(X,3);
    % show figure;
    if showEach~=0, figure; end
    for i=1:nodeNum
        mat = [];
        for j=1:gNum
            mat = [mat, X(i,:,j).'];
        end
        if showEach~=0
            figure;
        else
            subplot(nodeNum,1,i);
        end
        plot(mat); ylim(yRange);
        if i~=nodeNum, xticklabels({}); end
    end
end
