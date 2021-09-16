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
    c = (X.'*Y);
    if any(X) && any(Y)
        s = c/sqrt((X.'*X)*(Y.'*Y));
    else
        s = 0;
    end
end
