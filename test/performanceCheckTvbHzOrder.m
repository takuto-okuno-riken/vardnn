function performanceCheckTvbHzOrder
    node_num = 8;
    hzs = [64, 128, 256, 512, 1024, 2048];
    scans = [20, 20, 16, 8, 4, 2];
    pat = 1;
    orders = [2, 4, 8, 16, 32, 64];
    Gth = 0.2;

    h = length(hzs);
    od = length(orders);
    
    fcAUC = nan(h,od);
    gcAUC = nan(h,od);
    linueAUC = nan(h,od);

    fcRf = figure;
    %dlRf = figure;

    for i=1:length(hzs)
        hz = hzs(i);
        num_scan = scans(i);
        tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(pat) '-' num2str(hz) 'hz.mat'];
        load(tvbFile);

        % show original signal FC
        FC = calcFunctionalConnectivity(si);
        figure(fcRf); hold on; [~,~, fcAUC(i,1)] = plotROCcurve(FC, weights, 100, 1, Gth); hold off;
        title('ROC curve of FC');
        
        gcRf = figure;
        linueRf = figure;
        for j=1:length(orders)
            lag = orders(j);
            % show original signal granger causality index (mvGC)
            gcI = calcMultivariateGCI(si,lag);
            figure(gcRf); hold on; [~,~, gcAUC(i,j)] = plotROCcurve(gcI, weights, 100, 1, Gth); hold off;
            title(['ROC curve of GC (tvb sampling ' num2str(hz) 'hz)']);
            
            % linue TE result
            linueFile = ['results/tvb-hzodr-linue/linue_MultivAnalysis_tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(pat) '-' num2str(hz) 'hz-' num2str(lag) '.mat'];
            load(linueFile);
            A = outputToStore.reshapedMtx.';

            % show ROC curve of TE(LIN UE)
            figure(linueRf); hold on; [~,~, linueAUC(i,j)] = plotROCcurve(A, weights, 100, 1, Gth); hold off;        
            title(['ROC curve of LINUE-TE (tvb sampling ' num2str(hz) 'hz)']);
        end
    end
    figure; clims = [0.5, 1]; imagesc(fcAUC,clims); colorbar; daspect([1 1 1]);
    title('AUC map of FC');
    figure; clims = [0.5, 1]; imagesc(gcAUC,clims); colorbar; daspect([1 1 1]);
    title('AUC map of GC');
    figure; clims = [0.5, 1]; imagesc(linueAUC,clims); colorbar; daspect([1 1 1]);
    title('AUC map of LINUE-TE');

    fname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(pat) '-HzOrder_result.mat'];
    save(fname, 'fcAUC', 'gcAUC', 'linueAUC');
end

