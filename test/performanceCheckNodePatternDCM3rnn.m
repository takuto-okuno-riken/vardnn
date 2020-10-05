% Before using this function, download SPM12 codes from
% https://www.fil.ion.ucl.ac.uk/spm/software/download/
% and add a path "spm12" and sub folders, then remove "spm12/external" folder and sub folders.

% this script should run after performanceCheckNodePatternDCM3, performanceCheckNodePatternDCM3d and RNN-GC result
function performanceCheckNodePatternDCM3rnn
    % -------------------------------------------------------------------------
    N  = 8;
    T  = 300;                             % number of observations (scans)
    n  = 8;                               % number of regions or nodes

    prefix = 'net-pat3-';                 % original weight file prefix (result of *NodePatternDCM3d.m)
    Gth = 0;                              % 0 for pat3. 0.2 for pat4.
%    prefix = 'net-pat4-';                  % original weight file prefix (result of *NodePatternDCM3d.m)
%    Gth = 0.2;                             % 0 for pat3. 0.2 for pat4.

    %% pattern 1 -------------------------------------------------
%%{
    disp('network density 0.05');
    checkingPattern(N,T,n,prefix,Gth,1);
%%}
    %% pattern 2 -------------------------------------------------
%%{
    disp('network density 0.11');
    checkingPattern(N,T,n,prefix,Gth,2);
%%}
    %% pattern 6 -------------------------------------------------
%%{
    disp('network density 0.16');
    checkingPattern(N,T,n,prefix,Gth,6);
%%}
    %% pattern 7 -------------------------------------------------
%%{
    disp('network density 0.27');
    checkingPattern(N,T,n,prefix,Gth,7);
%%}
    %% pattern 8 -------------------------------------------------
%%{
    disp('network density 0.36');
    checkingPattern(N,T,n,prefix,Gth,8);
%%}
end

%% 
function checkingPattern(N,T,n,prefix,Gth,idx)
    fname = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) 'result.mat'];
    load(fname);

    % init
    rnnAUC = zeros(1,N);
    linueAUC = zeros(1,N);
    nnnueAUC = zeros(1,N);
    rnnROC = cell(N,2);
    linueROC = cell(N,2);
    nnnueROC = cell(N,2);
    rnnRf = figure;
    linueRf = figure;
    nnnueRf = figure;
    origf = figure;
    rnnTrial = 8;

    % reading RNN-GC, TE(LIN UE), TE(BIN NUE) results
    for k=1:N
        dlcmFile = ['results/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(dlcmFile);

        % show original connection
        figure(origf); plotDcmEC(pP.A);

        % -----------------------------------------------------------------
        % read RNN-GC result
        A = zeros(n,n);
        for j=1:rnnTrial
            rnnFile = ['results/rnngc/' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' sprintf('%02d',k) '_' num2str(j) '.txt'];
            Aj = readmatrix(rnnFile);
            A = A + Aj.';
            %figure; plotDcmEC(Aj.');
        end
        A = A / rnnTrial;
        %figure; plotDcmEC(A);

        % show ROC curve of RNN-GC
        figure(rnnRf); hold on; [rnnROC{k,1}, rnnROC{k,2}, rnnAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;
        
        % -----------------------------------------------------------------
        % read Transfer Entropy (LIN UE) result
        linueFile = ['results/linue/linue_MultivAnalysis_' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(linueFile);
        A = outputToStore.reshapedMtx.';

        % show ROC curve of TE(LIN UE)
        figure(linueRf); hold on; [linueROC{k,1}, linueROC{k,2}, linueAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;        

        % -----------------------------------------------------------------
        % read Transfer Entropy (NearestNeighber NUE) result
        nnnueFile = ['results/nnnue/nnnue_MultivAnalysis_' prefix num2str(n) 'x' num2str(T) '-idx' num2str(idx) '-' num2str(k) '.mat'];
        load(nnnueFile);
        A = outputToStore.reshapedMtx.';

        % show ROC curve of TE(NearestNeighber NUE)
        figure(nnnueRf); hold on; [nnnueROC{k,1}, nnnueROC{k,2}, nnnueAUC(k)] = plotROCcurve(A, pP.A, 100, 1, Gth); hold off;        
    end
    % save result
    save(fname, 'fcAUC','pcAUC','wcsAUC','gcAUC','dlAUC','dlgAUC','dcmAUC','rnnAUC','linueAUC','nnnueAUC', 'fcROC','pcROC','wcsROC','gcROC','dlROC','dlgROC','dcmROC','rnnROC','linueROC','nnnueROC');
    
    % show average ROC curve of DCM
    figure; 
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(wcsROC, N, [0.9,0.5,0]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
    plotErrorROCcurve(dcmROC, N, [0.2,0.2,0.8]);
    plotErrorROCcurve(rnnROC, N, [0.8,0.8,0.2]);
    plotErrorROCcurve(linueROC, N, [0.2,0.6,0.8]);
    plotErrorROCcurve(nnnueROC, N, [0.8,0.2,0.8]);
    plotErrorROCcurve(dlgROC, N, [0.6,0.6,0.6]);
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '--', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(dlROC, N, '-', [0.2,0.2,0.2],1.2);
    plotAverageROCcurve(dcmROC, N, '-', [0.2,0.2,0.8],0.5);
    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5);
    plotAverageROCcurve(linueROC, N, '--', [0.2,0.5,0.7],0.5);
    plotAverageROCcurve(nnnueROC, N, '--', [0.7,0.2,0.7],0.5);
    plotAverageROCcurve(dlgROC, N, '--', [0.6,0.6,0.6],0.5);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve idx' num2str(idx)]);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

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
