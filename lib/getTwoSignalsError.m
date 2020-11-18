%%
% Get error between two signals
% input:
%  X          multivariate time series matrix (node x time series)
%  Y          multivariate time series matrix (node x time series)

function [mae, maeerr, errs] = getTwoSignalsError(X, Y)
    nodeNum = size(X,1);
    d=0; e=0; errs=[];
    for i=1:nodeNum
        Xi = X(i,:);
        Yi = Y(i,:);

        d = d + sqrt(immse(Yi, Xi));
        e = e + mean(abs(Yi - Xi));
        A = single(Yi - Xi);
        errs = [errs; A];
    end
    mae = e / nodeNum;
    maeerr = std(errs(:),1) / sqrt(length(errs(:))); % uncorrected
end
