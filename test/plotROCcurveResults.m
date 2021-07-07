% this script should run after performanceCheckNodePatternDCM3, performanceCheckNodePatternDCM3d, performanceCheckNodePatternDCM3rnn and RNN-GC result
function plotROCcurveResults
    % -------------------------------------------------------------------------
    N  = 8;
    T  = 300;                             % number of observations (scans)
    n  = 8;                               % number of regions or nodes

    prefix = 'net-pat3-';
    checkingPattern(N,T,n,prefix,1);
    checkingPattern(N,T,n,prefix,2);
    checkingPattern(N,T,n,prefix,6);
    checkingPattern(N,T,n,prefix,7);
    checkingPattern(N,T,n,prefix,8);

    prefix = 'net-pat4-';
    checkingPattern(N,T,n,prefix,1);
    checkingPattern(N,T,n,prefix,2);
    checkingPattern(N,T,n,prefix,6);
    checkingPattern(N,T,n,prefix,7);
    checkingPattern(N,T,n,prefix,8);
end

%% 
function checkingPattern(N,T,n,prefix,idx)
    fname = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    load(fname);
    
    % show comparison ROC curves result
    figure; 
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.5,0.1,0.1]);
    plotErrorROCcurve(pcpcROC, N, [0.5,0.1,0.1]);
    plotErrorROCcurve(lsopcROC, N, [0.5,0.1,0.1]); % SPC-EN
    plotErrorROCcurve(plspcROC, N, [0.5,0.1,0.1]);
    plotErrorROCcurve(pgcROC, N, [0.0,0.5,0.0]);
    plotErrorROCcurve(gcROC, N, [0.1,0.8,0.1]);
    plotErrorROCcurve(mpcvargcROC, N, [0.1,0.8,0.1]); % PCA-cGCM
    plotErrorROCcurve(mlsogcROC, N, [0.1,0.8,0.1]);
    plotErrorROCcurve(mplsgcROC, N, [0.1,0.8,0.1]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]); % VARDNN-GC
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % VARDNN-DI
    plotErrorROCcurve(pcgcROC, N, [0.7,0.7,0.2]); % PC-GC
    plotErrorROCcurve(rnnROC, N, [0.7,0.7,0.2]); % RNN-GC
    plotErrorROCcurve(linueROC, N, [0.2,0.5,0.8]); % LINUE-TE
    plotErrorROCcurve(nnnueROC, N, [0.2,0.5,0.8]); % NNNUE-TE
    
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '-', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pcpcROC, N, '--', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(lsopcROC, N, '-.', [0.5,0.1,0.1],0.5); % SPC-EN
    plotAverageROCcurve(plspcROC, N, ':', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '-', [0.0,0.5,0.0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(mpcvargcROC, N, '--', [0.1,0.8,0.1],0.5); % PCA-cGCM
    plotAverageROCcurve(mlsogcROC, N, '-.', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(mplsgcROC, N, ':', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5); % RNN-GC
    plotAverageROCcurve(pcgcROC, N, '-.', [0.7,0.7,0.2],0.5); % PC-GC
    plotAverageROCcurve(linueROC, N, '-', [0.2,0.5,0.8],0.5); % LINUE-TE
    plotAverageROCcurve(nnnueROC, N, '--', [0.2,0.5,0.8],0.5); % NNNUE-TE
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

