%%
% Prot two signals
% input:
%  X          multivariate time series matrix (node x time series)
%  Y          multivariate time series matrix (node x time series)

function [R] = plotTwoSignalsCorrelation(X, Y)
    s = scatter(X(:),Y(:),3,'filled');
    s.MarkerEdgeColor = [0.65,0.65,0.65];
    s.MarkerFaceColor = [0.65,0.65,0.65];
    daspect([1 1 1]);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    R = corr2(X(:),Y(:));
end
