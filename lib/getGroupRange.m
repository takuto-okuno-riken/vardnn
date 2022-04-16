%%
% get range value [min, max, mean, std] of whole group
% input:
%  CX              cells of multivariate time series matrix {node x time series}
function [gMin, gMax, gM, gStd] = getGroupRange(CX)
    A = [];
    for i=1:length(CX)
        A = [A, CX{i}];
    end
    gMin = min(A(:));
    gMax = max(A(:));
    gM = nanmean(A(:));
    gStd = nanstd(A(:),1);
end
