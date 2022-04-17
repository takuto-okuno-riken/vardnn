%%
% get range struct [min, max, mean, std] of whole group
% input:
%  CX              cells of multivariate time series matrix {node x time series}
function gRange = getGroupRange(CX)
    A = [];
    for i=1:length(CX)
        A = [A, CX{i}];
    end
    gRange.min = min(A(:));
    gRange.max = max(A(:));
    gRange.m = nanmean(A(:));
    gRange.s = nanstd(A(:),1);
end
