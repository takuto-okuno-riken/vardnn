%%
% plot DLCM network weight by each node
% input:
%  netDLCM   trained DLCM network
%  type      graph type ('bar' (default), 'whisker')

function plotDlcmWeight(netDLCM, type)
    if nargin < 2, type = 'bar'; end
    % weight plot
    disp('plot DLCM connection weight amang nodes');
    nodeNum = length(netDLCM.nodeLayers);
    for i=1:nodeNum
        weight = netDLCM.nodeNetwork{i, 1}.Layers(2, 1).Weights;
        mweight = mean(weight,1);
        eweight = std(weight,1) / sqrt(size(weight,1));
        switch type
            case 'whisker'
                % box-and-whisker plot
                figure;
                boxplot(weight);
            otherwise
                % bar plot
                figure;
                bar(mweight);
                hold on
                er = errorbar(1:length(mweight),mweight,eweight);    
                er.Color = [0 0 0];                            
                er.LineStyle = 'none';  
                hold off
        end
    end
end
