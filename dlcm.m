%%
% DLCM command line tool

function dlcm(varargin)

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
    handles.inControls = {};
    handles.commandError = 0;
    handles.dlgc = 0;
    handles.mvgc = 0;
    handles.pwgc = 0;
    handles.fc = 0;
    handles.te = 0;
    handles.lag = 3;
    handles.pval = 0;
    handles.transform = 0;
    handles.alpha = 0.05;
    handles.showSig = 0;
    handles.showMat = 0;
    handles.maxEpochs = 1000;
    handles.L2Regularization = 0.05;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-d','--dlgc'}
                handles.dlgc = 1;
            case {'-m','--mvgc'}
                handles.mvgc = 1;
            case {'-p','--pwgc'}
                handles.pwgc = 1;
            case {'-f','--fc'}
                handles.fc = 1;
            case {'-t','--te'}
                handles.te = 1;
            case {'--lag'}
                handles.lag = str2num(varargin{i+1});
                i = i + 1;
            case {'--pval'}
                handles.pval = 1;
            case {'--transform'}
                handles.transform = str2num(varargin{i+1});
                i = i + 1;
            case {'--epoch'}
                handles.maxEpochs = str2num(varargin{i+1});
                i = i + 1;
            case {'--l2'}
                handles.L2Regularization = str2num(varargin{i+1});
                i = i + 1;
            case {'--ex'}
                handles.exoFiles = split(varargin{i+1},':');
                i = i + 1;
            case {'--nctrl'}
                handles.nodeControls = split(varargin{i+1},':');
                i = i + 1;
            case {'--ectrl'}
                handles.inControls = split(varargin{i+1},':');
                i = i + 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showmat'}
                handles.showMat = 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'-v','--version'}
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
    disp('  -d, --dlgc          output DLCM Granger Causality matrix result (<filename>_dlgc.csv)');
    disp('  -m, --mvgc          output multivaliate Granger Causality matrix result (<filename>_mvgc.csv)');
    disp('  -p, --pwgc          output pair-wised Granger Causality matrix result (<filename>_pwgc.csv)');
    disp('  -t, --te            output (LINUE) Transfer Entropy matrix result (<filename>_te.csv)');
    disp('  -f, --fc            output Functional Conectivity matrix result (<filename>_fc.csv)');
    disp('  --lag num           time lag <num> for mvGC, pwGC and TE (default:3)');
    disp('  --pval              output p-value matrix of DLCM-GC, mvGC, pwGC, TE and FC (<filename>_*_pval.csv)');
    disp('  --transform type    signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --ex files          exogenouse input signal <files> for DLCM-GC (file1.csv[:file2.csv:...])');
    disp('  --nctrl files       node status control <files> for DLCM-GC (file1.csv[:file2.csv:...])');
    disp('  --ectrl files       exogenous input control <files> for DLCM-GC (file1.csv[:file2.csv:...])');
    disp('  --epoch num         training epoch number <num> for DLCM-GC (default:1000)');
    disp('  --l2 num            training L2Regularization <num> for DLCM-GC (default:0.05)');
    disp('  --showsig           show status signals of <filename>.csv');
    disp('  --showmat           show result matrix of DLCM-GC, mvGC, pwGC, TE and FC');
    disp('  -v, --version       show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    for i = 1:length(handles.csvFiles)
        % load status signal csv file
        fname = handles.csvFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
        T = readtable(fname);
        X = table2array(T);
        nodeNum = size(X,1);
        sigLen = size(X,2);
        [exedir,name,ext] = fileparts(fname);

        % load exogenous signal csv file
        inNum = 0;
        inSignal = [];
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
            inSignal = table2array(T);
            inNum = size(inSignal, 1);
            if size(inSignal,2) < sigLen
                disp(['error : exogenous signal length is smaller than node status signal length : ' exoname]);
                continue;
            end
            inSignal = inSignal(:,1:sigLen);
        end

        % load node control csv file
        nodeControl = [];
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
        if inNum > 0
            inControl = eye(nodeNum, inNum);
        else
            inControl = [];
        end
        if ~isempty(handles.inControls)
            if length(handles.inControls)==1
                incntrolname = handles.inControls{1};
            elseif length(handles.inControls) >= i
                incntrolname = handles.inControls{i};
            else
                disp(['error : bad exogenous input control file list with ' fname]);
                continue;
            end
            T = readtable(ndcntrolname);
            inControl = table2array(T);
        end

        % signal transform raw or not
        if handles.transform == 1
            [X, sig, m, maxsi, minsi] = convert2SigmoidSignal(X);
            if ~isempty(inSignal)
                [inSignal, sig, m, maxsi, minsi] = convert2SigmoidSignal(inSignal);
            end
        end

        % show inputing status signals
        if handles.showSig > 0
            figure; plot(X.');
            title(['Status Signals : ' name]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
        
        % calc DLCM Granger causality
        if handles.dlgc > 0
            % layer parameters
            netDLCM = initDlcmNetwork(X, inSignal, nodeControl, inControl);
            % training DLCM network
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
            netDLCM = trainDlcmNetwork(X, inSignal, nodeControl, inControl, netDLCM, options);
            [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
            disp(['end training : rsme=' num2str(rsme)]);
            
            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [gcI, h, P, F, cvFd, nodeAIC, nodeBIC] = plotDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM, 0, 0, handles.alpha)
                title(['DLCM Granger Causality Index : ' name]);
            else
                [gcI, h, P, F, cvFd, nodeAIC, nodeBIC] = calcDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM, handles.alpha);
            end
            % output result matrix csv file
            outputCsvFile(gcI, [name '_dlgc.csv']);

            % output result p-value matrix csv file
            if handles.pval > 0
                outputCsvFile(P, [name '_dlgc_pval.csv']);
            end
        end
        
        % calc multivaliate Granger causality
        if handles.mvgc > 0
            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [gcI, h, P, F, cvFd, nodeAIC, nodeBIC] = plotMultivariateGCI2(X, handles.lag, 0, 0, handles.alpha);
                title(['multivariate Granger Causality Index : ' name]);
            else
                [gcI, h, P, F, cvFd, nodeAIC, nodeBIC] = calcMultivariateGCI2(X, handles.lag, handles.alpha);
            end
            % output result matrix csv file
            outputCsvFile(gcI, [name '_mvgc.csv']);

            % output result p-value matrix csv file
            if handles.pval > 0
                outputCsvFile(P, [name '_mvgc_pval.csv']);
            end
        end
        
        % calc pair-wised Granger causality
        if handles.pwgc > 0
            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [gcI, h, P, F, cvFd, AIC, BIC] = plotPairwiseGCI(X, handles.lag, 0, 0, handles.alpha);
                title(['pairwised Granger Causality Index : ' name]);
            else
                [gcI, h, P, F, cvFd, AIC, BIC] = calcPairwiseGCI(X, handles.lag, handles.alpha);
            end
            % output result matrix csv file
            outputCsvFile(gcI, [name '_pwgc.csv']);

            % output result p-value matrix csv file
            if handles.pval > 0
                outputCsvFile(P, [name '_pwgc_pval.csv']);
            end
        end
        
        % calc Transfer Entropy (LINUE)
        if handles.te > 0
            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [TE, h, P, F, cvFd, nodeAIC, nodeBIC] = plotLinueTE(X, handles.lag, 0, 0, handles.alpha);
                title(['Transfer Entropy (LINUE) : ' name]);
            else
                [TE, h, P, F, cvFd, nodeAIC, nodeBIC] = calcLinueTE(X, handles.lag, handles.alpha);
            end
            % output result matrix csv file
            outputCsvFile(TE, [name '_te.csv']);

            % output result p-value matrix csv file
            if handles.pval > 0
                outputCsvFile(P, [name '_te_pval.csv']);
            end
        end
        
        % calc Function connectivity
        if handles.fc > 0
            % show original signal FC
            if handles.showMat > 0
                figure; [FC,P] = plotFunctionalConnectivity(X);
                title(['Functional Connectivity : ' name]);
            else
                [FC,P] = calcFunctionalConnectivity(X);
            end
            % output result matrix csv file
            outputCsvFile(FC, [name '_fc.csv']);
            
            % output result p-value matrix csv file
            if handles.pval > 0
                outputCsvFile(P, [name '_fc_pval.csv']);
            end
        end
    end
end

%%
% output csv file function
function outputCsvFile(mat, outfname)
    T = array2table(mat);
    writetable(T,outfname,'WriteVariableNames',false);
    disp(['output csv file : ' outfname]);
end
