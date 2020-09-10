function performanceCheckNodePatternTVBovfnn
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
    trial = 6;

    % init
    dlAUC = zeros(N,trial);
    dlErr = zeros(N,trial);
    origf = figure;
    origSigf = figure;

    for k=1:N
        tvbFile = ['data/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-' num2str(k) '.mat'];
        load(tvbFile);
        density = length(find(weights>Gth)) / (node_num * (node_num-1));

        [si, sig, m, maxsi, minsi] = convert2SigmoidSignal(si);
        [uu, sig2, m2, maxsi2, minsi2] = convert2SigmoidSignal(uu);
            
        % show original connection
        figure(origf); plotDcmEC(weights);
        figure(origSigf); plot(t, si);

        % calcurate and show DLCM-GC
        nodeNum = size(si,1);
        sigLen = size(si,2);
        inputNum = size(uu,1);
        inControl = eye(nodeNum, nodeNum);
        for j=1:trial
            dlcmFile = ['results/net-patrww-'  num2str(nodeNum) 'x' num2str(num_scan) '-idx' num2str(i) '-' num2str(k) 'ovfnn' num2str(j) '.mat'];
            if exist(dlcmFile, 'file')
                load(dlcmFile);
            else
                % train DLCM
                Y = si;
                inSignal = uu;

                % estimate neuron number of hidden layers
                hiddenNums = zeros(2,1);
                hiddenNums(1) = 10 * j;
                hiddenNums(2) = ceil(hiddenNums(1)*2/3);

                % set initial bias for each neuron
                biasMat = ones(hiddenNums(1),1) * 0;

                % layer parameters
                netDLCM = createDlcmNetwork(nodeNum, inputNum, hiddenNums, nodeControl, inControl, [], [], biasMat);

                % training DLCM network
                maxEpochs = 1000;
                miniBatchSize = ceil(sigLen / 3);
                options = trainingOptions('adam', ...
                    'ExecutionEnvironment','cpu', ...
                    'MaxEpochs',maxEpochs, ...
                    'MiniBatchSize',miniBatchSize, ...
                    'Shuffle','every-epoch', ...
                    'GradientThreshold',5,...
                    'L2Regularization',0.05, ...
                    'Verbose',false);
            %            'Plots','training-progress');

                disp('start training');
                netDLCM = trainDlcmNetwork(Y, inSignal, [], inControl, netDLCM, options);
                save(dlcmFile, 'netDLCM', 'Y', 'inSignal', 'Y', 'sig', 'm', 'maxsi', 'minsi', 'sig2', 'm2', 'maxsi2', 'minsi2');
            end
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);
            dlErr(k,j) = rsme;
            
            % show DLCM-GC
            dlGC = calcDlcmGCI(Y, inSignal, [], inControl, netDLCM);

            % calc ROC curve
            [~, ~, dlAUC(k,j)] = calcROCcurve(dlGC, weights, 100, 1, Gth);
        end
    end

    % save result
    fname = ['results/tvb-wongwang' num2str(node_num) 'x' num2str(num_scan) 'scan-pat' num2str(i) '-' num2str(hz) 'hz-ovfnn-result.mat'];
    save(fname, 'dlAUC','dlErr');
end
