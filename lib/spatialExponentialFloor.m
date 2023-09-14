%%
% calculate spatialExponentialFloor from tools (Shinn et al., 2022)
% returns SE (node x node)
% input:
%  D           distance matrix (node x node)
%  saLambda    SA-λ value of spatial autocorrelation
%  saInf       SA-∞ value of spatial autocorrelation

function SE = spatialExponentialFloor(D, saLambda, saInf)
    SE = exp(-D/saLambda)*(1-saInf)+saInf;
end
