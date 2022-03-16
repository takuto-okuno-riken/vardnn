%%
% VARDNN command line tool

function vardnn(varargin)

    % set version number
    versionNumber = '0.1';

    % add script path
    if ~isdeployed % checking MATLAB mode or stand-alone mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind).file;
        [exedir,exename,ext] = fileparts(relpath);
        if exist([exedir '/util'],'dir')
            addpath([exedir '/util']);
            addpath([exedir '/lib']);
        end
    end

    % get exe file full path
    global exePath;
    global exeName;
    [exePath, exeName, ext] = exeFilename();

    % init command line input
    handles.commandError = 0;
    handles.csvFiles = {};
    handles.exoFiles = {};
    handles.nodeControls = {};
    handles.exControls = {};
    handles.groundTruth = {};
    handles.roiNames = {};
    handles.commandError = 0;
    handles.outpath = 'results';
    handles.vddi = 0;
    handles.vdgc = 0;
    handles.mvgc = 0;
    handles.pwgc = 0;
    handles.te = 0;
    handles.fc = 0;
    handles.pc = 0;
    handles.wc = 0;
    handles.mvar = 0;
    handles.dllag = 1;
    handles.lag = 3;
    handles.pval = 0;
    handles.fval = 0;
    handles.aic = 0;
    handles.bic = 0;
    handles.transform = 0;
    handles.transopt = NaN;
    handles.format = 0;
    handles.alpha = 0.05;
    handles.Gth = 0;
    handles.showSig = 0;
    handles.showEx = 0;
    handles.showMat = 0;
    handles.showCG = 0;
    handles.showROC = 0;
    handles.maxEpochs = 1000;
    handles.L2Regularization = 0.05;
    handles.noCache = 0;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-d','--vddi'}
                handles.vddi = 1;
            case {'-c','--vdgc'}
                handles.vdgc = 1;
            case {'-m','--mvgc'}
                handles.mvgc = 1;
            case {'-g','--pwgc'}
                handles.pwgc = 1;
            case {'-t','--te'}
                handles.te = 1;
            case {'-f','--fc'}
                handles.fc = 1;
            case {'-p','--pc'}
                handles.pc = 1;
            case {'-w','--wc'}
                handles.wc = 1;
            case {'-v','--var'}
                handles.mvar = 1;
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--lag'}
                handles.lag = str2num(varargin{i+1});
                i = i + 1;
            case {'--pval'}
                handles.pval = 1;
            case {'--fval'}
                handles.fval = str2num(varargin{i+1});
                i = i + 1;
            case {'--aic'}
                handles.aic = 1;
            case {'--bic'}
                handles.bic = 1;
            case {'--transform'}
                handles.transform = str2num(varargin{i+1});
                i = i + 1;
            case {'--transopt'}
                handles.transopt = str2num(varargin{i+1});
                i = i + 1;
            case {'--format'}
                handles.format = str2num(varargin{i+1});
                i = i + 1;
            case {'--epoch'}
                handles.maxEpochs = str2num(varargin{i+1});
                i = i + 1;
            case {'--l2'}
                handles.L2Regularization = str2num(varargin{i+1});
                i = i + 1;
            case {'--ex'}
                handles.exoFiles = strsplit(varargin{i+1},':');
                i = i + 1;
            case {'--nctrl'}
                handles.nodeControls = strsplit(varargin{i+1},':');
                i = i + 1;
            case {'--ectrl'}
                handles.exControls = strsplit(varargin{i+1},':');
                i = i + 1;
            case {'--groundtruth'}
                handles.groundTruth = strsplit(varargin{i+1},':');
                i = i + 1;
            case {'--roiname'}
                handles.roiNames = strsplit(varargin{i+1},':');
                i = i + 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showex'}
                handles.showEx = 1;
            case {'--showmat'}
                handles.showMat = 1;
            case {'--showcg'}
                handles.showCG = 1;
            case {'--showroc'}
                handles.showROC = 1;
            case {'--nocache'}
                handles.noCache = 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'--version'}
                disp([exeName ' version : ' num2str(versionNumber)]);
                return;
            otherwise
                if strcmp(varargin{i}(1), '-')
                    disp(['bad option : ' varargin{i}]);
                    i = size(varargin, 2);
                    handles.commandError = 1;
                else
                    handles.csvFiles = [handles.csvFiles varargin{i}];
                end
        end
        i = i + 1;
    end
    
    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif isempty(handles.csvFiles)
        disp('no input files. please specify node status signal files.');
        showUsage();
        return;
    end

    % process input files
    processInputFiles(handles);
end

%%
% show usage function
function showUsage()
    global exePath;
    global exeName;
    disp(['usage: ' exeName ' [options] filename.csv ...']);
    disp('  -d, --vddi          output VARDNN Directional Influence matrix result (<filename>_vddi.csv)');
    disp('  -c, --vdgc          output VARDNN Granger Causality matrix result (<filename>_vdgc.csv)');
    disp('  -m, --mvgc          output multivaliate Granger Causality matrix result (<filename>_mvgc.csv)');
    disp('  -g, --pwgc          output pairwise Granger Causality matrix result (<filename>_pwgc.csv)');
    disp('  -t, --te            output (LINUE) Transfer Entropy matrix result (<filename>_te.csv)');
    disp('  -f, --fc            output Functional Conectivity matrix result (<filename>_fc.csv)');
    disp('  -p, --pc            output Partial Correlation matrix result (<filename>_pc.csv)');
    disp('  -w, --wc            output Wavelet Coherence matrix result (<filename>_wc.csv)');
    disp('  -v, --var           output VAR (Vector Auto-Regression) coefficient matrix result (<filename>_mvar.csv)');
    disp('  --outpath           output files path (default:"results")');
    disp('  --pval              save P-value matrix of VARDNN-GC, mvGC, pwGC, TE, FC and PC (<filename>_*_pval.csv)');
    disp('  --fval alpha        save F-value with <alpha> matrix of VARDNN-GC, mvGC, pwGC and TE (<filename>_*_fval.csv, <filename>_*_fcrit.csv)');
    disp('  --aic               save AIC matrix of VARDNN-GC, mvGC, pwGC and TE (<filename>_*_aic.csv)');
    disp('  --bic               save BIC matrix of VARDNN-GC, mvGC, pwGC and TE (<filename>_*_bic.csv)');
    disp('  --format type       save file format <type> 0:csv, 1:mat(each), 2:mat(all) (default:0)');
    disp('  --groundtruth files calculate ROC curve and save AUC of VARDNN-DI, VARDNN-GC, mVAR, mvGC, pwGC, TE, FC, PC and WC (<filename>_*_auc.csv)');
    disp('  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --transopt num      signal transform option <num> (for type 1:centroid value)');
    disp('  --lag num           time lag <num> for mvGC, pwGC, TE and mVAR (default:3)');
    disp('  --ex files          VARDNN exogenous input signal <files> (file1.csv[:file2.csv:...])');
    disp('  --nctrl files       VARDNN node status control <files> (file1.csv[:file2.csv:...])');
    disp('  --ectrl files       VARDNN exogenous input control <files> (file1.csv[:file2.csv:...])');
    disp('  --epoch num         VARDNN training epoch number <num> (default:1000)');
    disp('  --l2 num            VARDNN training L2Regularization <num> (default:0.05)');
    disp('  --roiname files     ROI names <files> (file1.csv[:file2.csv:...])');
    disp('  --showsig           show node status signals of <filename>.csv');
    disp('  --showex            show exogenous input signals of <file1>.csv');
    disp('  --showmat           show result matrix of VARDNN-DI, VARDNN-GC, mVAR, mvGC, pwGC, TE, FC, PC and WC');
    disp('  --showcg            show circle graph of VARDNN-DI, VARDNN-GC, mVAR, mvGC, pwGC, TE, FC, PC and WC');
    disp('  --showroc           show ROC curve (by GroundTruth) of VARDNN-DI, VARDNN-GC, mVAR, mvGC, pwGC, TE, FC, PC and WC');
    disp('  --nocache           do not use cache file for VARDNN training');
    disp('  --version           show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;

    % init
    N = length(handles.csvFiles);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    wcROC = cell(N,2);
    gcROC = cell(N,2);
    pwROC = cell(N,2);
    dlROC = cell(N,2);
    dleROC = cell(N,2);
    teROC = cell(N,2);
    mvarROC = cell(N,2);
    if handles.showROC > 0
        if handles.vddi > 0, dleRf = figure; end
        if handles.vdgc > 0, dlRf = figure; end
        if handles.mvgc > 0, gcRf = figure; end
        if handles.pwgc > 0, pwRf = figure; end
        if handles.te > 0, teRf = figure; end
        if handles.fc > 0, fcRf = figure; end
        if handles.pc > 0, pcRf = figure; end
        if handles.wc > 0, wcRf = figure; end
        if handles.mvar > 0, mvarecRf = figure; end
    end
    
    % process each file
    for i = 1:N
        % init data
        X = [];
        exSignal = [];
        nodeControl = [];
        exControl = [];
        inNum = 0;
        groundTruth = [];
        roiNames = [];
        auc = NaN;

        % load node status signals csv or mat file
        fname = handles.csvFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
        [path,name,ext] = fileparts(fname);
        if strcmp(ext,'.mat')
            load(fname);
            inNum = size(exSignal, 1);
        else
            T = readtable(fname);
            X = table2array(T);
        end
        nodeNum = size(X,1);
        sigLen = size(X,2);
        
        if handles.format==2 % if save format is mat(all)
            if i==1
                savename = name;
            end
        else
            savename = name;
        end

        % load exogenous input signals csv file
        if ~isempty(handles.exoFiles)
            if length(handles.exoFiles)==1
                exoname = handles.exoFiles{1};
            elseif length(handles.exoFiles) >= i
                exoname = handles.exoFiles{i};
            else
                disp(['error : bad exogenous file list with ' fname]);
                continue;
            end
            T = readtable(exoname);
            exSignal = table2array(T);
            inNum = size(exSignal, 1);
            if size(exSignal,2) < sigLen
                disp(['error : exogenous signal length is smaller than node status signal length : ' exoname]);
                continue;
            end
            exSignal = exSignal(:,1:sigLen);
        end

        % load node control csv file
        if ~isempty(handles.nodeControls)
            if length(handles.nodeControls)==1
                ndcntrolname = handles.nodeControls{1};
            elseif length(handles.nodeControls) >= i
                ndcntrolname = handles.nodeControls{i};
            else
                disp(['error : bad node control file list with ' fname]);
                continue;
            end
            T = readtable(ndcntrolname);
            nodeControl = table2array(T);
        end

        % load exogenous input control csv file
        if inNum > 0 && isempty(exControl)
            exControl = ones(nodeNum, inNum);
        end
        if ~isempty(handles.exControls)
            if length(handles.exControls)==1
                incntrolname = handles.exControls{1};
            elseif length(handles.exControls) >= i
                incntrolname = handles.exControls{i};
            else
                disp(['error : bad exogenous input control file list with ' fname]);
                continue;
            end
            T = readtable(incntrolname);
            exControl = table2array(T);
        end

        % load ground truth csv file
        if ~isempty(handles.groundTruth)
            if length(handles.groundTruth)==1
                gtname = handles.groundTruth{1};
            elseif length(handles.groundTruth) >= i
                gtname = handles.groundTruth{i};
            else
                disp(['error : bad ground truth file list with ' fname]);
                continue;
            end
            T = readtable(gtname);
            groundTruth = table2array(T);
        end

        % load ROI names csv file
        roiNames = generateNumberROINames(nodeNum+inNum);
        if ~isempty(handles.roiNames)
            if length(handles.roiNames)==1
                roiname = handles.roiNames{1};
            elseif length(handles.roiNames) >= i
                roiname = handles.roiNames{i};
            else
                disp(['error : bad ROI name file list with ' fname]);
                continue;
            end
            T = readtable(roiname, 'ReadVariableNames', 0);
            tn = table2array(T);
            roiNames(1:length(tn)) = tn;
        end

        % signal transform raw or not
        if handles.transform == 1
            [X, sig, c, maxsi, minsi] = convert2SigmoidSignal(X, handles.transopt);
            if ~isempty(exSignal)
                [exSignal, sig, c, maxsi, minsi] = convert2SigmoidSignal(exSignal, handles.transopt);
            end
        end

        % show node status signals
        if handles.showSig > 0
            figure; plot(X.');
            title(['Node Status Signals : ' name]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
        % show exogenous input signals
        if handles.showEx > 0 && ~isempty(exSignal)
            figure; plot(exSignal.');
            title(['Exogenous Input Signals : ' exoname]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
        
        % train VARDNN network
        if handles.vddi > 0 || handles.vdgc > 0
            vardnnFile = [exePath '/results/cache-vardnn-' name '.mat'];
            if exist(vardnnFile, 'file') && handles.noCache == 0
                disp(['read cache file : ' vardnnFile]);
                load(vardnnFile);
            else
                % layer parameters
                net = initMvarDnnNetwork(X, exSignal, nodeControl, exControl, handles.dllag);
                % training VARDNN network
                miniBatchSize = ceil(sigLen / 3);
                options = trainingOptions('adam', ...
                    'ExecutionEnvironment','cpu', ...
                    'MaxEpochs', handles.maxEpochs, ...
                    'MiniBatchSize', miniBatchSize, ...
                    'Shuffle', 'every-epoch', ...
                    'GradientThreshold', 5,...
                    'L2Regularization', handles.L2Regularization, ...
                    'Verbose',false);

                disp('start training');
                net = trainMvarDnnNetwork(X, exSignal, nodeControl, exControl, net, options);
                [time, loss, rsme] = getMvarDnnTrainingResult(net);
                disp(['VARDNN training result : rsme=' num2str(rsme)]);
                if handles.noCache == 0
                    save(vardnnFile, 'net');
                end
            end
        end
        
        % groundTruth size check
        isFullNode = 1;
        if ~isempty(groundTruth) && size(groundTruth,2) == nodeNum
            isFullNode = 0;
        end

        % calc VARDNN-DI
        if handles.vddi > 0
            % show VARDNN-DI matrix
            if handles.showMat > 0
                figure; DI = plotMvarDnnDI(net, nodeControl, exControl, 0, isFullNode);
                title(['VARDNN(' num2str(handles.dllag) ') Directional Influence : ' name]);
            else
                DI = calcMvarDnnDI(net, nodeControl, exControl, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(DI, ['VARDNN(' num2str(handles.dllag) ') Directional Influence : ' name], roiNames, 1, 5, 1);
            end
            
            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(dleRf); hold on; [dleROC{i,1}, dleROC{i,2}, auc] = plotROCcurve(DI, groundTruth, 100, 1, handles.Gth); hold off;
                    title(['ROC curve of VARDNN(' num2str(handles.dllag) ') Directional Influence']);
                else
                    [dleROC{i,1}, dleROC{i,2}, auc] = calcROCcurve(DI, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, DI, [], [], [], [], [], auc, [savename '_vddi']);
        end

        % calc VARDNN-GC
        if handles.vdgc > 0
            % show VARDNN-GC matrix
            if handles.showMat > 0
                figure; [dlGC, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotMvarDnnGCI(X, exSignal, nodeControl, exControl, net, 0, handles.alpha, isFullNode);
                title(['VARDNN(' num2str(handles.dllag) ') Granger Causality Index : ' name]);
            else
                [dlGC, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMvarDnnGCI(X, exSignal, nodeControl, exControl, net, handles.alpha, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(dlGC, ['VARDNN(' num2str(handles.dllag) ') Granger Causality Index : ' name], roiNames, 1, 5, 1);
            end
            
            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(dlRf); hold on; [dlROC{i,1}, dlROC{i,2}, auc] = plotROCcurve(dlGC, groundTruth, 100, 1, handles.Gth); hold off;
                    title(['ROC curve of VARDNN(' num2str(handles.dllag) ') Granger Causality']);
                else
                    [dlROC{i,1}, dlROC{i,2}, auc] = calcROCcurve(dlGC, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, dlGC, P, F, cvFd, AIC, BIC, auc, [savename '_vdgc']);
        end
        
        % calc multivaliate Granger causality
        if handles.mvgc > 0
            % show mvGC matrix
            if handles.showMat > 0
                figure; [mvGC, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotMultivariateGCI(X, exSignal, nodeControl, exControl, handles.lag, 0, handles.alpha, isFullNode);
                title(['multivariate Granger Causality(' num2str(handles.lag) ') Index : ' name]);
            else
%                mvGCa = calcMultivariateGCI_(X, exSignal, nodeControl, exControl, handles.lag, isFullNode);
%                mvGCb = calcMultivariateGCI__(X, exSignal, nodeControl, exControl, handles.lag, isFullNode);
                [mvGC, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMultivariateGCI(X, exSignal, nodeControl, exControl, handles.lag, handles.alpha, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(mvGC, ['multivariate Granger Causality(' num2str(handles.lag) ') Index : ' name], roiNames, 1, 5, 1);
            end
            
            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(gcRf); hold on; [gcROC{i,1}, gcROC{i,2}, auc] = plotROCcurve(mvGC, groundTruth, 100, 1, handles.Gth); hold off;
                    title(['ROC curve of multivariate Granger Causality(' num2str(handles.lag) ')']);
                else
                    [gcROC{i,1}, gcROC{i,2}, auc] = calcROCcurve(mvGC, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, mvGC, P, F, cvFd, AIC, BIC, auc, [savename '_mvgc']);
        end
        
        % calc pair-wised Granger causality
        if handles.pwgc > 0
            % show pwGC matrix
            if handles.showMat > 0
                figure; [pwGC, h, P, F, cvFd, AIC, BIC] = plotPairwiseGCI(X, exSignal, nodeControl, exControl, handles.lag, 0, handles.alpha, isFullNode);
                title(['pairwised Granger Causality(' num2str(handles.lag) ') Index : ' name]);
            else
                [pwGC, h, P, F, cvFd, AIC, BIC] = calcPairwiseGCI(X, exSignal, nodeControl, exControl, handles.lag, handles.alpha, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(pwGC, ['pairwised Granger Causality(' num2str(handles.lag) ') Index : ' name], roiNames, 1, 5, 1);
            end

            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(pwRf); hold on; [pwROC{i,1}, pwROC{i,2}, auc] = plotROCcurve(pwGC, groundTruth, 100, 1, handles.Gth); hold off;
                    title(['ROC curve of pairwised Granger Causality(' num2str(handles.lag) ')']);
                else
                    [pwROC{i,1}, pwROC{i,2}, auc] = calcROCcurve(pwGC, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, pwGC, P, F, cvFd, AIC, BIC, auc, [savename '_pwgc']);
        end
        
        % calc Transfer Entropy (LINUE)
        if handles.te > 0
            % show TE matrix
            if handles.showMat > 0
                figure; [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotLinueTE(X, exSignal, nodeControl, exControl, handles.lag, 0, handles.alpha, isFullNode);
                title(['Transfer Entropy(' num2str(handles.lag) ') (LINUE) : ' name]);
            else
                [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcLinueTE(X, exSignal, nodeControl, exControl, handles.lag, handles.alpha, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(TE, ['Transfer Entropy(' num2str(handles.lag) ') (LINUE) : ' name], roiNames, 1, 5, 1);
            end

            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(teRf); hold on; [teROC{i,1}, teROC{i,2}, auc] = plotROCcurve(TE, groundTruth, 100, 1, handles.Gth); hold off;
                    title(['ROC curve of Transfer Entropy(' num2str(handles.lag) ') (LINUE)']);
                else
                    [teROC{i,1}, teROC{i,2}, auc] = calcROCcurve(TE, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, TE, P, F, cvFd, AIC, BIC, auc, [savename '_te']);
        end
        
        % calc Function connectivity
        if handles.fc > 0
            % show FC matrix
            if handles.showMat > 0
                figure; [FC,P] = plotFunctionalConnectivity(X, exSignal, nodeControl, exControl, isFullNode);
                title(['Functional Connectivity : ' name]);
            else
                [FC,P] = calcFunctionalConnectivity(X, exSignal, nodeControl, exControl, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(FC, ['Functional Connectivity : ' name], roiNames, 0.4, 0.8, 0.2, true);
            end

            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(fcRf); hold on; [fcROC{i,1}, fcROC{i,2}, auc] = plotROCcurve(FC, groundTruth, 100, 1, handles.Gth); hold off;
                    title('ROC curve of Functional Connectivity');
                else
                    [fcROC{i,1}, fcROC{i,2}, auc] = calcROCcurve(FC, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, FC, P, [], [], [], [], auc, [savename '_fc']);
        end

        % calc Partial Correlation
        if handles.pc > 0
            % show PC matrix
            if handles.showMat > 0
                figure; [PC,P] = plotPartialCorrelation(X, exSignal, nodeControl, exControl, isFullNode);
                title(['Partial Correlation : ' name]);
            else
                [PC,P] = calcPartialCorrelation(X, exSignal, nodeControl, exControl, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(PC, ['Partial Correlation : ' name], roiNames, 0.4, 0.8, 0.2, true);
            end
            
            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(pcRf); hold on; [pcROC{i,1}, pcROC{i,2}, auc] = plotROCcurve(PC, groundTruth, 100, 1, handles.Gth); hold off;
                    title('ROC curve of Partial Correlation');
                else
                    [pcROC{i,1}, pcROC{i,2}, auc] = calcROCcurve(PC, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, PC, P, [], [], [], [], auc, [savename '_pc']);
        end

        % calc Wavelet Coherence
        if handles.wc > 0
            % show WC matrix
            if handles.showMat > 0
                figure; [mWCS, WCOH, WCS] = plotWaveletCoherence(X, exSignal, nodeControl, exControl, isFullNode);
                title(['Wavelet Coherence : ' name]);
            else
                [mWCS, WCOH, WCS] = calcWaveletCoherence(X, exSignal, nodeControl, exControl, isFullNode);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(mWCS, ['Wavelet Coherence : ' name], roiNames, 0.4, 0.8, 0.2, true);
            end

            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(wcRf); hold on; [wcROC{i,1}, wcROC{i,2}, auc] = plotROCcurve(mWCS, groundTruth, 100, 1, handles.Gth); hold off;
                    title('ROC curve of Wavelet Coherence');
                else
                    [wcROC{i,1}, wcROC{i,2}, auc] = calcROCcurve(mWCS, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, mWCS, [], [], [], [], [], auc, [savename '_wc']);
        end

        % calc mVAR
        if handles.mvar > 0
            netMVAR = initMvarNetwork(X, exSignal, nodeControl, exControl, handles.lag);
            % show MVAR matrix
            [~, ~, mvar] = calcMvarDI(netMVAR, nodeControl, exControl, isFullNode);
            if handles.showMat > 0
                figure; plotDirectedFC(mvar, 'mVAR', 5);
                title(['mVAR(' num2str(handles.lag) ') coefficient : ' name]);
            end
            
            if handles.showCG > 0
                figure; plotCircleGraph(mvar, ['mVAR(' num2str(handles.lag) ') coefficient : ' name], roiNames, 1, 5, 1);
            end
            
            % show ROC curve 
            if ~isempty(groundTruth)
                if handles.showROC > 0
                    figure(mvarecRf); hold on; [mvarROC{i,1}, mvarROC{i,2}, auc] = plotROCcurve(mvar, groundTruth, 100, 1, handles.Gth); hold off;
                    title(['ROC curve of mVAR(' num2str(handles.lag) ') coefficient']);
                else
                    [mvarROC{i,1}, mvarROC{i,2}, auc] = calcROCcurve(mvar, groundTruth, 100, 1, handles.Gth);
                end
            end
            
            % output result matrix files
            saveResultFiles(handles, mvar, [], [], [], [], [], auc, [savename '_mvar']);
        end
    end
    
    % show average ROC curve
    if handles.showROC > 0
        figure; 
        hold on;
        plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
        plotErrorROCcurve(pcROC, N, [0.8,0.2,0.2]);
        plotErrorROCcurve(wcROC, N, [0.8,0.5,0.2]);
        plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
        plotErrorROCcurve(pwROC, N, [0.0,0.5,0.0]);
        plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
        plotErrorROCcurve(dleROC, N, [0.2,0.2,0.2]);
        plotErrorROCcurve(teROC, N, [0.2,0.6,0.8]);
        plotErrorROCcurve(mvarROC, N, [0.3,0.3,0.3]);
        plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
        plotAverageROCcurve(pcROC, N, '--', [0.8,0.2,0.2],0.5);
        plotAverageROCcurve(wcROC, N, '--', [0.8,0.5,0.2],0.5);
        plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
        plotAverageROCcurve(pwROC, N, '--', [0.0,0.5,0.0],0.5);
        plotAverageROCcurve(dlROC, N, '--', [0.2,0.2,0.2],0.8);
        plotAverageROCcurve(dleROC, N, '-', [0.2,0.2,0.2],1.2);
        plotAverageROCcurve(teROC, N, '--', [0.2,0.5,0.7],0.5);
        plotAverageROCcurve(mvarROC, N, '-', [0.3,0.3,0.3],0.5);
        plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
        hold off;
        ylim([0 1]);
        xlim([0 1]);
        daspect([1 1 1]);
        title(['averaged ROC curve of ' num2str(N) ' data']);
        xlabel('False Positive Rate')
        ylabel('True Positive Rate')
    end
end

%%
function roiNames = generateNumberROINames(nodeNum)
    roiNames = cell(1, nodeNum);
    for i=1:nodeNum
        roiNames{i} = num2str(i);
    end
end

%%
% show cercle graph
function plotCircleGraph(mat, matname, roiNames, rmin, rmax, rstep, isRaw)
    if nargin < 7
        isRaw = 0;
    end
    if isRaw
        mOrg = mat;
    else
        sigma = std(mat(:),1,'omitnan');
        avg = mean(mat(:),'omitnan');
        mOrg = (mat - avg) / sigma;
    end
    rangeW = -rmin:-rstep:-rmax;
    rangeS = rmin:rstep:rmax;
    
    % fix matrix size to Square
    if size(mOrg,1) < size(mOrg,2)
        m = nan(size(mOrg,2));
        m(1:size(mOrg,1),:) = mOrg;
        mOrg = m;
    end
    % plot weaker weights
    cstep = 0.8 / length(rangeW);
    for i=1:length(rangeW)
        m = mOrg;
        m(m>rangeW(i)) = 0;
        if i<length(rangeW)
            m(m<=rangeW(i+1)) = 0;
        end
        m(isnan(m)) = 0;
        hold on; 
        G = digraph(m', 'omitselfloops');
        gp=plot(G);
        layout(gp,'circle');
        gp.LineStyle = '-';
        gp.NodeLabel = [];
        gp.EdgeColor = [1-i*cstep, 1-i*cstep, 0.9];
        hold off;
    end
    % plot stronger weights
    cstep = 0.8 / length(rangeS);
    for i=1:length(rangeS)
        m = mOrg;
        m(m<rangeS(i)) = 0;
        if i<length(rangeS)
            m(m>=rangeS(i+1)) = 0;
        end
        m(isnan(m)) = 0;
        hold on; 
        G = digraph(m', 'omitselfloops');
        gp=plot(G);
        layout(gp,'circle');
        gp.LineStyle = '-';
        gp.NodeLabel = [];
        gp.EdgeColor = [0.9, 1-i*cstep, 1-i*cstep];
        hold off;
    end
    gp.NodeColor = [0.7, 0.7, 0.7];
    gp.NodeLabel = roiNames(1:size(mOrg,1));
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    daspect([1 1 1]);
    title(['Cricle Graph of ' matname]);
end

%%
% output result matrix files
%
function saveResultFiles(handles, Index, P, F, cvFd, AIC, BIC, auc, outname)
    if handles.format == 1
        save([handles.outpath '/' outname '.mat'],'Index', 'P', 'F', 'cvFd', 'AIC', 'BIC', 'auc');
    elseif handles.format == 2
        fname = [handles.outpath '/' outname '_all.mat'];
        if exist(fname,'file')
            t = load(fname);
            t.Index(:,:,end+1) = Index; Index = t.Index;
            t.P(:,:,end+1) = P; P = t.P;
            t.F(:,:,end+1) = F; F = t.F;
            t.cvFd(:,:,end+1) = cvFd; cvFd = t.cvFd;
            t.AIC(:,:,end+1) = AIC; AIC = t.AIC;
            t.BIC(:,:,end+1) = BIC; BIC = t.BIC;
            t.auc(end+1) = auc; auc = t.auc;
        end
        save(fname,'Index', 'P', 'F', 'cvFd', 'AIC', 'BIC', 'auc');
    else
        % output result matrix csv file
        outputCsvFile(Index, [handles.outpath '/' outname '.csv']);

        % output result P-value matrix csv file
        if handles.pval > 0 && ~isempty(P)
            outputCsvFile(P, [handles.outpath '/' outname '_pval.csv']);
        end

        % output result F-value matrix csv file
        if handles.fval > 0 && ~isempty(F)
            outputCsvFile(F, [handles.outpath '/' outname '_fval.csv']);
            outputCsvFile(cvFd, [handles.outpath '/' outname '_fcrit.csv']);
        end

        % output AIC matrix csv file
        if handles.aic > 0 && ~isempty(AIC)
            outputCsvFile(AIC, [handles.outpath '/' outname '_aic.csv']);
        end

        % output BIC matrix csv file
        if handles.bic > 0 && ~isempty(BIC)
            outputCsvFile(BIC, [handles.outpath '/' outname '_bic.csv']);
        end

        % output auc csv file
        if ~isnan(auc)
            outputCsvFile(auc, [handles.outpath '/' outname '_auc.csv']);
        end
    end
end

%%
% output csv file function
%
function outputCsvFile(mat, outfname)
    T = array2table(mat);
    writetable(T,outfname,'WriteVariableNames',false);
    disp(['output csv file : ' outfname]);
end
