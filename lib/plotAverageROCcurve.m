
function [x, y] = plotAverageROCcurve(roc, N, line, col, width)
    maxCol = 0;
    for k=1:N
        if size(roc{k,1},2)>maxCol, maxCol=size(roc{k,1},2); end
    end
    x = nan(N,maxCol);
    y = nan(N,maxCol);
    for k=1:N
        len = size(roc{k,1},2);
        if len==0, continue; end
        x(k,:) = roc{k,1}(end);
        y(k,:) = roc{k,2}(end);
        x(k,1:len) = roc{k,1}(1:len);
        y(k,1:len) = roc{k,2}(1:len);
    end
    mx = mean(x,1);
    my = mean(y,1);
    % plot line
    plot(mx, my, line,'Color',col,'LineWidth',width);
end
