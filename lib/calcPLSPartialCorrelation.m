%%
% Caluclate PLS Partial Correlation
% returns PLS Partial Correlation (PC)
% input:
%  X       multivariate time series matrix (node x time series)
%  lambda  regularisation level

function [PC] = calcPLSPartialCorrelation(X)
    n = size(X,1);
    PC = nan(n,n);
    ncomp = floor(n / 5);
    if ncomp < 3, ncomp = 3; end
    if ncomp > 50, ncomp = 50; end
    Y = 1:n;
    for i=1:n
        for j=i:n
            x = X(i,:).';
            y = X(j,:).';
            Yij = setdiff(setdiff(Y,i),j);
            z = X(Yij,:).';

            [XL,YL,XS,YS,b,PCTVAR,MSE,stats1] = plsregress(z,x,ncomp);
            [XL,YL,XS,YS,b,PCTVAR,MSE,stats2] = plsregress(z,y,ncomp);
            r1 = stats1.Yresiduals;
            r2 = stats2.Yresiduals;
            PC(i,j) = (r1.'*r2) / (sqrt(r1.'*r1)*sqrt(r2.'*r2));
            PC(j,i) = PC(i,j);
        end
    end
end
