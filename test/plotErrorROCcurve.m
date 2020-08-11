
function [x, y] = plotErrorROCcurve(roc, N, col)
    x = zeros(N,size(roc{1,1},2));
    y = zeros(N,size(roc{1,1},2));
    for k=1:N
        x(k,:) = roc{k,1};
        y(k,:) = roc{k,2};
    end
    mx = mean(x,1);
    my = mean(y,1);
    % error area
    errx = std(x,0,1) / sqrt(size(x,1));
    erry = std(y,0,1) / sqrt(size(y,1));
    xt = mx -errx;
    yt = my +erry;
    xb = mx +errx;
    yb = my -erry;
    Y = [yb fliplr(yt)];
    X = [xb fliplr(xt)];
    patch(X,Y,col,'EdgeColor','none','FaceAlpha',.04);
end