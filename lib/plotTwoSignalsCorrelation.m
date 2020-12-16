%%
% Prot two signals
% input:
%  X          multivariate time series matrix (node x time series)
%  Y          multivariate time series matrix (node x time series)
%  color      edge and face color
%  marker     marker type
%  sz         marker size

function [R] = plotTwoSignalsCorrelation(X, Y, color, marker, sz)
    if nargin < 5
        sz = 3;
    end
    if nargin < 4
        marker = 'o';
    end
    if nargin < 3
        color = [0.65,0.65,0.65];
    end
    s = scatter(X(:),Y(:),sz,'filled',marker);
    s.MarkerEdgeColor = color;
    s.MarkerFaceColor = color;
    daspect([1 1 1]);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    R = corr2(X(:),Y(:));
end
