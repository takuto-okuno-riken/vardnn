%%
% Caluclate PartiallyConditionedGrangerCausality
% returns estimated causality (GC)
% input:
%  X            multivariate time series matrix (node x time series)
%  lags         number of lags for autoregression (default:1)
%  ndRate       (of ndmax) as the value selected as the knee of the curve (default:0.8)

% Before using this function, download PartiallyConditionedGrangerCausality codes from
% https://github.com/danielemarinazzo/PartiallyConditionedGrangerCausality
% and add a path "PartiallyConditionedGrangerCausality-master" and sub folders. 

function GC = calcPCGC(X, lags, ndRate)
    if nargin < 3, ndRate = 0.8; end
    if nargin < 2, lags = 1; end
    nvar = size(X,1);
    
    % ndmax calculation
    % ndmax , the maximum number of variables to consider as simultaneously
    % contributing. As a rule of thumb you can choose ndmax = nvar/2 for nvar<100
    % and reduce this fraction to nvar/20 for 1000 regions, but this of course depends
    % on the nature of the system.
    if nvar < 20, ndmax = nvar-1;
    elseif nvar < 100, ndmax = floor(nvar/(1+nvar*0.01));
    else ndmax = floor(nvar/(2+nvar*0.02));
    end
    
    [y, ind]=init_partial_conditioning_par_m(X.',ndmax,lags);
%    plot(1:nodeNum-2,diff(y'));
    
    pcgc=partial_CGC_fix_nd_m(X.',lags,ceil(ndmax*ndRate),ind);
    GC = pcgc.';
end
