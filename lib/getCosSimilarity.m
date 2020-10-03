%%
% Get cosign between two signals
% input:
%  X          any matrix
%  Y          any matrix

function s = getCosSimilarity(X, Y)
    X2 = X(:);
    Y2 = Y(:);
    X = X2(~isnan(X2) & ~isnan(Y2));
    Y = Y2(~isnan(Y2) & ~isnan(X2));
    s = (X.'*Y)/sqrt((X.'*X)*(Y.'*Y));
end
