%%
% Get cosign between two signals
% input:
%  X          any matrix
%  Y          any matrix

function s = getCosSimilarity(X, Y)
    X = X(:);
    Y = Y(:);
    s = (X.'*Y)/sqrt((X.'*X)*(Y.'*Y));
end
