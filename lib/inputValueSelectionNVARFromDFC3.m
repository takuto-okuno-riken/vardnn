%%
% Get control matrix from input directed FC
% input:
%   dFC       3D directed FC matrix (node x node x lags)
%   rate      top n% for active connection [0,1]

function control3 = inputValueSelectionNVARFromDFC3(dFC, rate)
    control3 = zeros(1,size(dFC,2),size(dFC,3));
    dFC2 = mean(dFC,1);
    [~,idx] = sort(abs(dFC2(:)),'descend');
    imax = ceil(length(idx) * rate);
    control3(idx(1:imax)) = 1;
end
