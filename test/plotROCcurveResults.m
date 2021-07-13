% this script should run after performanceCheckNodePatternDCM3, performanceCheckNodePatternDCM3d, performanceCheckNodePatternDCM3rnn and RNN-GC result
function plotROCcurveResults
    % -------------------------------------------------------------------------
    N  = 8;
    T  = 300;                             % number of observations (scans)
    n  = 8;                               % number of regions or nodes

    idxs = [1,2,6,7,8];
    for i=1:length(idxs)
        idx = idxs(i);
        fname = ['results/net-pat3-' num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
        plotROCresuts(N,fname,idx);
    end

    for i=1:length(idxs)
        idx = idxs(i);
        fname = ['results/net-pat4-' num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
        plotROCresuts(N,fname,idx);
    end
    
    node_nums = [16,32,48,64,80,98];
    num_scan = 55;
    hz = 64;
    N = 8;
    for i=1:length(node_nums)
        fname = ['results/tvb-wongwang' num2str(node_nums(i)) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-result.mat'];
        plotROCresuts(N,fname,i);
    end
    
    % AD vs HC
    fname = ['results/ad/ad-cn-ad-roi' num2str(132) '-result.mat'];
    plotROCresuts(60,fname,1);
end

%% 
function plotROCresuts(N,fname,idx)
    load(fname);
    
    if exist('teAUC','var'), linueAUC=teAUC; linueROC=teROC; end
    X = [fcAUC.', pcAUC.', pcpcAUC.', lsopcAUC.', plspcAUC.', pgcAUC.', gcAUC.', mpcvargcAUC.', mlsogcAUC.', mplsgcAUC.'];
    if exist('rnnROC','var'), X = [X, rnnAUC.']; else X = [X, zeros(length(fcAUC),1)]; end
    X = [X, pcgcAUC.', linueAUC.'];
    if exist('nnnueROC','var'), X = [X, nnnueAUC.']; else X = [X, zeros(length(fcAUC),1)]; end
    X = [X, dlAUC.', dlwAUC.', pcdlAUC.', pcdlwAUC.'];
    
    % show comparison ROC curves result
    figure; 
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.5,0.1,0.1]);
    plotErrorROCcurve(pcpcROC, N, [0.5,0.1,0.1]);
    plotErrorROCcurve(lsopcROC, N, [0.5,0.1,0.1]); % SPC-EN
    plotErrorROCcurve(plspcROC, N, [0.5,0.1,0.1]);
    plotErrorROCcurve(pgcROC, N, [0.2,0.4,0.2]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(mpcvargcROC, N, [0.2,0.8,0.2]); % PCA-cGCM
    plotErrorROCcurve(mlsogcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(mplsgcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]); % VARDNN-GC
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % VARDNN-DI
    plotErrorROCcurve(pcdlROC, N, [0.5,0.5,0.5]); % PCVARDNN-GC
    plotErrorROCcurve(pcdlwROC, N, [0.5,0.5,0.5]); % PCVARDNN-DI
    plotErrorROCcurve(pcgcROC, N, [0.7,0.7,0.2]); % PC-GC
    if exist('rnnROC','var'), plotErrorROCcurve(rnnROC, N, [0.7,0.7,0.2]); end % RNN-GC
    plotErrorROCcurve(linueROC, N, [0.2,0.5,0.8]); % LINUE-TE
    if exist('nnnueROC','var'), plotErrorROCcurve(nnnueROC, N, [0.2,0.5,0.8]); end % NNNUE-TE
    
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '-', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pcpcROC, N, '--', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(lsopcROC, N, '-.', [0.5,0.1,0.1],0.5); % SPC-EN
    plotAverageROCcurve(plspcROC, N, ':', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '-', [0.2,0.4,0.2],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.2,0.8,0.2],0.5);
    plotAverageROCcurve(mpcvargcROC, N, '--', [0.2,0.8,0.2],0.5); % PCA-cGCM
    plotAverageROCcurve(mlsogcROC, N, '-.', [0.2,0.8,0.2],0.5);
    plotAverageROCcurve(mplsgcROC, N, ':', [0.2,0.8,0.2],0.5);
    if exist('rnnROC','var'), plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5); end % RNN-GC
    plotAverageROCcurve(pcgcROC, N, '-.', [0.7,0.7,0.2],0.5); % PC-GC
    plotAverageROCcurve(linueROC, N, '-', [0.2,0.5,0.8],0.5); % LINUE-TE
    if exist('nnnueROC','var'), plotAverageROCcurve(nnnueROC, N, '--', [0.2,0.5,0.8],0.5); end % NNNUE-TE
    plotAverageROCcurve(pcdlROC, N, '--', [0.5,0.5,0.5],0.8); % PCVARDNN-GC
    plotAverageROCcurve(pcdlwROC, N, '-', [0.5,0.5,0.5],0.8); % PCVARDNN-DI
    plotAverageROCcurve(dlROC, N, '--', [0.2,0.2,0.2],0.8); % VARDNN-GC
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.2); % VARDNN-DI
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(idx)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end
