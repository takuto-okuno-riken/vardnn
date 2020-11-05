% load healthy control analysis result of DLCM-GC
load('results/ad-dlcm-cn-roi132.mat');
load('data/roiNames.mat');

% invert weight 
m = 1 ./ meanWeights;
G = digraph(m, 'omitselfloops');

% plot graph
gp=plot(G);
layout(gp,'force','WeightEffect','direct');
gp.NodeLabel = roiNames;
gp.LineStyle = ':';
gp.NodeColor = [1 0 0];
gp.EdgeColor = [0.9 0.9 0.9];


% plot graph
sigma = std(meanWeights(:),'omitnan');
avg = mean(meanWeights(:),'omitnan');
mOrg = (meanWeights - avg) / sigma;
rangeW = [-2,-3,-4,-5];
rangeS = [2,3,4,5];
figure;
for i=1:length(rangeW)
    m = mOrg;
    m(m>rangeW(i)) = 0;
    if i<length(rangeW)
        m(m<=rangeW(i+1)) = 0;
    end
    hold on; 
    G = digraph(m, 'omitselfloops');
    gp=plot(G);
    layout(gp,'circle');
    gp.LineStyle = '-';
    gp.EdgeColor = [1-i*0.2, 1-i*0.2, 0.9];
    hold off;
end
for i=1:length(rangeS)
    m = mOrg;
    m(m<rangeS(i)) = 0;
    if i<length(rangeS)
        m(m>=rangeS(i+1)) = 0;
    end
    hold on; 
    G = digraph(m, 'omitselfloops');
    gp=plot(G);
    layout(gp,'circle');
    gp.LineStyle = '-';
    gp.EdgeColor = [0.9, 1-i*0.2, 1-i*0.2];
    hold off;
end
gp.NodeColor = [0.7, 0.7, 0.7];
gp.NodeLabel = roiNames;
