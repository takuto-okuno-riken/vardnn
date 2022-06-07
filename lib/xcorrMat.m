% XCORR Cross-correlation function estimates (two matrix version)
% x and y must have same column number
% based on matlab xcorr()

function [c,lags] = xcorrMat(x,y,maxlag,scale)

    % Perform the cross-correlation.
    c = crosscorr(x,y,maxlag);
    
    % Scale the output.
    c = scaleOutput(scale,c,x,y);
    
    if nargout == 2
        lags = -maxlag:maxlag;
    end
end

%--------------------------------------------------------------------------
% Compute cross-correlation for matrix inputs.
% x and y must have same column number
function c = crosscorr(x,y,maxlag)
    [m,n] = size(x);
    mxl = min(maxlag,m - 1);

    m2 = findTransformLength(m);
    X = fft(x,m2,1);
    Y = fft(y,m2,1);
    C = reshape(X,m2,1,n).*conj(Y(:,:));
    if isreal(x) && isreal(y)
        c1 = ifft(C,[],1,'symmetric');
    else
        c1 = ifft(C,[],1);
    end
    % c1 is M2-by-N-by-N.
    % Keep only the lags we want, and move the negative lags before the
    % positive lags. Also flatten the result to 2-D.
    c = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];
end

%--------------------------------------------------------------------------
function m = findTransformLength(m)
    m = 2*m;
    while true
        r = m;
        for p = [2 3 5 7]
            while (r > 1) && (mod(r, p) == 0)
                r = r / p;
            end
        end
        if r == 1
            break;
        end
        m = m + 1;
    end
end

%--------------------------------------------------------------------------
% Scale correlation as specified.
function c = scaleOutput(scale,c,x,y)
    if strcmp(scale,'none')
        return
    end

    m = size(x,1);
    if strcmp(scale,'biased')
        % Scales the raw cross-correlation by 1/M.
        c = c./m;
    elseif strcmp(scale,'unbiased')
        % Scales the raw correlation by 1/(M-abs(lags)).
        L = (size(c,1) - 1)/2;
        scaleUnbiased = (m - abs(-L:L)).';
        scaleUnbiased(scaleUnbiased <= 0) = 1;
        c = c./scaleUnbiased;
    else % 'normalized'/'coeff'
        % Compute autocorrelations at zero lag.
        % scale = norm(x)*norm(y) is numerically superior but slower.
        cxx0 = sqrt(sum(x.*x));
        cyy0 = sqrt(sum(y.*y));
        cxx0(cxx0==0) = 1e-15; % avoid divide by 0
        cyy0(cyy0==0) = 1e-15; % avoid divide by 0
        scaleCoeffCross = cyy0'*cxx0;
        c = c./scaleCoeffCross(:)';
    end
end
