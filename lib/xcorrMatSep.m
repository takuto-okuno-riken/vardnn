% XCORR Auto-correlation function estimates (matrix recursive version)
% x and y must have same column number
% based on matlab xcorr()

function [c,lags] = xcorrMatSep(x,y,maxlag,scale,sepNum)
    if sepNum <= 1
        [c,lags] = xcorrMat(x,y,maxlag,scale);
    else
        % separate matrix
        n = size(x,2);
        s = ceil(n / sepNum);
        cx = cell(sepNum,1);
        cy = cell(sepNum,1);
        for i=1:sepNum
            if i==sepNum
                cx{i} = x(:,1+s*(i-1):end);
                cy{i} = y(:,1+s*(i-1):end);
            else
                cx{i} = x(:,1+s*(i-1):s*i);
                cy{i} = y(:,1+s*(i-1):s*i);
            end
        end
        clear x; clear y;
        % calc xcorr
        CB = cell(sepNum,sepNum);
        for i=1:sepNum
            for j=i:sepNum
                [c,lags] = xcorrMat(cx{i},cy{j},maxlag,scale);
                B = zeros(size(cy{j},2),size(cx{i},2),length(lags));
                for k=1:length(lags)
                    B(:,:,k) = reshape(c(k,:),size(cy{j},2),size(cx{i},2));
                end
                CB{j,i} = B;
            end
        end
        clear cx; clear cy;

        % conjunct separate results
        B = zeros(n,n,length(lags));
        c = zeros(length(lags),n*n);
        sc = length(lags);
        for k=1:sc
            m = [];
            for i=1:sepNum
                b = [];
                for j=1:sepNum
                    if isempty(CB{j,i})
                        b = [b; CB{i,j}(:,:,sc+1-k)'];
                    else
                        b = [b; CB{j,i}(:,:,k)];
                    end
                end
                m = [m, b];
            end
            B(:,:,k) = m;
            c(k,:) = m(:);
        end
    end
end
