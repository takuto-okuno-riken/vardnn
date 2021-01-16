function performanceCheckNodePatternTVBovfl2
    num_scan = 55;
    if num_scan == 54  % oh's mouse 98 node. density around 0.15. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [1, 1, 1, 1, 1, 1];
    elseif num_scan == 55  % oh's mouse 98 node. density 0.15. weight add.
        node_nums = [16,32,48,64,80,98];
        Gths = [1, 1, 1, 1, 1, 1];
    end
    % test sparse and full density
    hz = 64;
    N = 8;

    for i=1:length(node_nums)
        checkingPattern(node_nums(i), num_scan, hz, Gths(i), N, i);
    end
end

function checkingPattern(node_num, num_scan, hz, Gth, N, i)
    l2rs = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
    trial = length(l2rs);

    % init
    dlAUC = zeros(N,trial);
    dlErr = zeros(N,trial);
    origf = figure;
    origSigf = figure;

    for k=1:N
        tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
        load(tvbFile);
        density = length(find(weights>Gth)) / (node_num * (node_num-1));

        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(si);
        [uu, sig2, c2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            
        % show original connection
        figure(origf); plotEC(weights, 'Ground Truth', 1);
        figure(origSigf); plot(t, si);

        % calcurate and show DLCM-GC
        nodeNum = size(si,1);
        sigLen = size(si,2);
        inControl = eye(nodeNum, nodeNum);
        for j=1:trial
            dlcmFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) 'ovfl' num2str(j) '.mat'];
            if exist(dlcmFile, 'file')
                load(dlcmFile);
            else
                % train DLCM
                Y = si;
                inSignal = uu;
                % layer parameters
                netDLCM = initDlcmNetwork(Y, inSignal, [], inControl);
                % training DLCM network
                maxEpochs = 1000;
                miniBatchSize = ceil(sigLen / 3);
                options = trainingOptions('adam', ...
                    'ExecutionEnvironment','cpu', ...
                    'MaxEpochs',maxEpochs, ...
                    'MiniBatchSize',miniBatchSize, ...
                    'Shuffle','every-epoch', ...
                    'GradientThreshold',5,...
                    'L2Regularization',l2rs(j), ...
                    'Verbose',false);
            %            'Plots','training-progress');

                disp('start training');
                netDLCM = trainDlcmNetwork(Y, inSignal, [], inControl, netDLCM, options);
                save(dlcmFile, 'netDLCM', 'Y', 'inSignal', 'Y', 'sig', 'c', 'maxsi', 'minsi', 'sig2', 'c2', 'maxsi2', 'minsi2');
            end
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);
            dlErr(k,j) = rsme;
            
            % show DLCM-GC
            dlGC = calcDlcmGCI_(Y, inSignal, [], inControl, netDLCM);

            % calc ROC curve
            [~, ~, dlAUC(k,j)] = calcROCcurve(dlGC, weights, 100, 1, Gth);
        end
    end

    % save result
    fname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-ovfl-result.mat'];
    save(fname, 'dlAUC','dlErr');
end
