%%
% Get control matrix from input directed FC
% input:
%   dFC       3D directed FC matrix (node x node x lags)
%   rate      top n% for active connection [0,1]

function control3 = inputValueSelectionFromDFC3(dFC, rate)
    control3 = zeros(size(dFC,1),size(dFC,2),size(dFC,3));
    [~,idx] = sort(abs(dFC(:)),'descend');
    imax = ceil(length(idx) * rate);
    control3(idx(1:imax)) = 1;
end
