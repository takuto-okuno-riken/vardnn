%%
% Calculate Z-scores of matrix (outmat)
% input:
%  mat    input matrix

function outmat = calcZScores(mat)
    sigma = nanstd(mat(:),1);
    avg = nanmean(mat(:));
    outmat = (mat - avg) / sigma;
end
