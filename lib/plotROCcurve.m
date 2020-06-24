%%
% Plot ROC curve of X with ground truth G
% returns x, y vectors and AUC value
% input:
%  X           target matrix (node x node) to get ROC curve
%  G           ground truth matrix (node x node) (TRUE is G > Gth)
%  ignoreDiag  ignore diagonal in the matrix (default:1)
%  Gth         ground truth threshold (default:0)

function [x, y, auc] = plotROCcurve(X, G, ignoreDiag, Gth)
    if nargin < 4, Gth = 0; end
    if nargin < 3, ignoreDiag = 1; end
    [x, y, auc] = calcROCcurve(X, G, ignoreDiag, Gth);
    plot(x, y);
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title('ROC curve');
end
