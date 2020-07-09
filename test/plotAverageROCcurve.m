
function [x, y] = plotAverageROCcurve(roc, N, line, col, width)
    x = zeros(N,size(roc{1,1},2));
    y = zeros(N,size(roc{1,1},2));
    for k=1:N
        x(k,:) = roc{k,1};
        y(k,:) = roc{k,2};
    end
    mx = mean(x,1);
    my = mean(y,1);
    % plot line
    plot(mx, my, line,'Color',col,'LineWidth',width);
end
