
function [mx, my] = plotAverageROCcurve(roc, N, line, col, width)
    step = 100;
    mx = 0:step;
    y = nan(N,length(mx));
    for k=1:N
        tX = roc{k,1};
        tY = roc{k,2};
        lasti = 1;
        for i = 1:length(tX)
            if i < length(tX) && tX(i) == tX(i+1), continue; end
            
            % find mx;
            mx1=floor(tX(lasti)*step);
            mx2=floor(tX(i)*step);
            idx1=find(mx==mx1);
            idx2=find(mx==mx2);
            % linear interpolation
            if tY(lasti) == tY(i) || idx1 == idx2
                y(k,idx1:idx2) = tY(i);
            else
                dy = tY(i)-tY(lasti);
                for j=idx1:idx2
                    y(k,j) = dy*((j-idx1)/(idx2-idx1)) + tY(lasti);
                end
            end
            lasti = i;
        end
        for j=2:length(mx)
            if isnan(y(k,j)), y(k,j)=y(k,j-1); end
        end
        y(k,end) = 1;
%        hold on; plot(mx, y(k,:)); hold off;
    end
    mx = mx / step;
    mx = [0, mx, 1];
    my = [0, mean(y,1), 1];
    % plot line
    plot(mx, my, line,'Color',col,'LineWidth',width);
end
