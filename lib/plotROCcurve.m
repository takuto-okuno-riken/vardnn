%%
% Plot ROC curve of X with ground truth G
% returns x, y vectors and AUC value
% input:
%  X           target matrix (node x node) to get ROC curve
%  G           ground truth matrix (node x node) (TRUE is G > Gth)
%  step        step number of ground truth (default:100)
%  ignoreDiag  ignore diagonal in the matrix (default:1)
%  Gth         ground truth threshold (default:0)

function [x, y, auc] = plotROCcurve(X, G, step, ignoreDiag, Gth)
    if nargin < 5, Gth = 0; end
    if nargin < 4, ignoreDiag = 1; end
    if nargin < 3, step = 100; end
    [x, y, auc] = calcROCcurve(X, G, step, ignoreDiag, Gth);
    hold on;
    plot(x, y);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;    
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title('ROC curve');
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
