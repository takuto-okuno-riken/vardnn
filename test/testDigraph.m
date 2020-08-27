% load healthy control analysis result of DLCM-GC
load('results/ad-dlcm-cn-roi132.mat');
load('data/roiNames.mat');

% invert weight 
m = 1 ./ meanWeights;
% G = digraph(m, 'upper', 'omitselfloops'); % for FC
G = digraph(m, 'omitselfloops');

% plot graph
gp=plot(G);
layout(gp,'force','WeightEffect','direct');
%layout(gp,'circle');
gp.NodeLabel = roiNames;
gp.LineStyle = ':';
gp.NodeColor = [1 0 0];
gp.EdgeColor = [0.9 0.9 0.9];
