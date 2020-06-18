%%
% Estimate initial weight for DLCM
% input:
%  sz      vector of (neuron num) (node num) 
%  scale   weight scale

function initWeight = estimateInitWeightRoughHe(sz, scale)
    filterSize = [sz(1) sz(2)];
    numIn = filterSize(1) * filterSize(2);
    varscale = 0.1;
    
    varWeights = 2 / ((1 + varscale^2) * numIn);
    initWeight = randn(sz) * sqrt(varWeights) * scale;
end