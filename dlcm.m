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
    handles.fval = 0;
    handles.aic = 0;
    handles.bic = 0;
    handles.transform = 0;
    handles.alpha = 0.05;
    handles.showSig = 0;
    handles.showEx = 0;
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
                handles.inControls = strsplit(varargin{i+1},':');
                i = i + 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showex'}
                handles.showEx = 1;
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
    disp('  --pval              output P-value matrix of DLCM-GC, mvGC, pwGC, TE and FC (<filename>_*_pval.csv)');
    disp('  --fval alpha        output F-value with <alpha> matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_fval.csv, <filename>_*_fcrit.csv)');
    disp('  --aic               output AIC matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_aic.csv)');
    disp('  --bic               output BIC matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_bic.csv)');
    disp('  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)');
    disp('  --lag num           time lag <num> for mvGC, pwGC and TE (default:3)');
    disp('  --ex files          DLCM exogenouse input signal <files> (file1.csv[:file2.csv:...])');
    disp('  --nctrl files       DLCM node status control <files> (file1.csv[:file2.csv:...])');
    disp('  --ectrl files       DLCM exogenous input control <files> (file1.csv[:file2.csv:...])');
    disp('  --epoch num         DLCM training epoch number <num> (default:1000)');
    disp('  --l2 num            DLCM training L2Regularization <num> (default:0.05)');
    disp('  --showsig           show node status signals of <filename>.csv');
    disp('  --showex            show exogenous input signals of <file1>.csv');
    disp('  --showmat           show result matrix of DLCM-GC, mvGC, pwGC, TE and FC');
    disp('  -v, --version       show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;
    for i = 1:length(handles.csvFiles)
        % load node status signals csv file
        fname = handles.csvFiles{i};
        if ~exist(fname,'file')
            disp(['file is not found. ignoring : ' fname]);
            continue;
        end
        T = readtable(fname);
        X = table2array(T);
        nodeNum = size(X,1);
        sigLen = size(X,2);
        [path,name,ext] = fileparts(fname);

        % load exogenous input signals csv file
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
            T = readtable(incntrolname);
            inControl = table2array(T);
        end

        % signal transform raw or not
        if handles.transform == 1
            [X, sig, m, maxsi, minsi] = convert2SigmoidSignal(X);
            if ~isempty(inSignal)
                [inSignal, sig, m, maxsi, minsi] = convert2SigmoidSignal(inSignal);
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
        if handles.showEx > 0 && ~isempty(inSignal)
            figure; plot(inSignal.');
            title(['Exogenous Input Signals : ' exoname]);
            xlabel('Time Series');
            ylabel('Signal Value');
        end
        
        % calc DLCM Granger causality
        if handles.dlgc > 0
            dlcmFile = [exePath '/results/cache-dlcm-' name '.mat'];
            if exist(dlcmFile, 'file')
                disp(['read cache file : ' dlcmFile]);
                load(dlcmFile);
            else
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
                disp(['DLCM training result : rsme=' num2str(rsme)]);
                save(dlcmFile, 'netDLCM');
            end

            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM, 0, 0, handles.alpha);
                title(['DLCM Granger Causality Index : ' name]);
            else
                [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcDlcmGCI(X, inSignal, nodeControl, inControl, netDLCM, handles.alpha);
            end
            
            % output result matrix files
            outputCsvFiles(handles, gcI, P, F, cvFd, AIC, BIC, [name '_dlgc']);
        end
        
        % calc multivaliate Granger causality
        if handles.mvgc > 0
            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotMultivariateGCI2(X, handles.lag, 0, 0, handles.alpha);
                title(['multivariate Granger Causality Index : ' name]);
            else
                [gcI, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcMultivariateGCI2(X, handles.lag, handles.alpha);
            end
            
            % output result matrix files
            outputCsvFiles(handles, gcI, P, F, cvFd, AIC, BIC, [name '_mvgc']);
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

            % output result matrix files
            outputCsvFiles(handles, gcI, P, F, cvFd, AIC, BIC, [name '_pwgc']);
        end
        
        % calc Transfer Entropy (LINUE)
        if handles.te > 0
            % show original signal granger causality index 
            if handles.showMat > 0
                figure; [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = plotLinueTE(X, handles.lag, 0, 0, handles.alpha);
                title(['Transfer Entropy (LINUE) : ' name]);
            else
                [TE, h, P, F, cvFd, AIC, BIC, nodeAIC, nodeBIC] = calcLinueTE(X, handles.lag, handles.alpha);
            end
            
            % output result matrix files
            outputCsvFiles(handles, TE, P, F, cvFd, AIC, BIC, [name '_te']);
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
            
            % output result matrix files
            outputCsvFiles(handles, FC, P, [], [], [], [], [name '_fc']);
        end
    end
end

%%
% output result matrix files
%
function outputCsvFiles(handles, mat, P, F, cvFd, AIC, BIC, outname)
    % output result matrix csv file
    outputCsvFile(mat, [outname '.csv']);

    % output result P-value matrix csv file
    if handles.pval > 0 && ~isempty(P)
        outputCsvFile(P, [outname '_pval.csv']);
    end

    % output result F-value matrix csv file
    if handles.fval > 0 && ~isempty(F)
        outputCsvFile(F, [outname '_fval.csv']);
        outputCsvFile(cvFd, [outname '_fcrit.csv']);
    end
    
    % output AIC matrix csv file
    if handles.aic > 0 && ~isempty(AIC)
        outputCsvFile(AIC, [outname '_aic.csv']);
    end
    
    % output BIC matrix csv file
    if handles.bic > 0 && ~isempty(BIC)
        outputCsvFile(BIC, [outname '_bic.csv']);
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
