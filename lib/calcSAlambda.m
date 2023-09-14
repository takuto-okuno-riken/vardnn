%%
% calculate the SA-λ and SA-∞ measures of spatial autocorrelation (Shinn et al., 2022)
% returns SA-λ, SA-∞
% input:
%  cm              NxN correlation matrix of timeseries (node x node)
%  dist            NxN distance matrix, representing the spatial distance (node x node) (output of pdist2 matrix)
%  discretization  The size of the bins to use when computing the SA parameters. (default: 1)
%                  Data that has values up to around 100 should be fine with the default. 

function [SAlambda, SAinf] = calcSAlambda(cm, dist, discretization)
    if nargin < 3, discretization = 1; end
    distFlat = round(dist(:) / discretization) * discretization;
    binX = unique(distFlat)';
    bl = length(binX);
    Y = discretize(distFlat(:),binX); % returns bl-1 bin group ids
    cmBin = nan(1,bl-1);
    for i=1:bl-1
        idx = find(Y==i);
        cmBin(i) = mean(cm(idx));
    end
%    binX(1) = 1; % ?? not sure
    binX = binX(1:end-1);
    v0 = [10, .3]; % initial value
    fun = @(v)sum((cmBin - spatialExponentialFloor(binX, v(1), v(2))).^2,'all'); % error evaluation function
%    options = optimset('Display','iter');
    [x,fval] = fminsearch(fun,v0); %,options);
    SAlambda = x(1); SAinf = x(2);
end
