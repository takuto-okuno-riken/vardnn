function analyzeAlzheimerDLCM
    % fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-70_F_CN_nii', 'ADNI2_65-70_M_CN_nii', 'ADNI2_70-75_F_CN_nii', 'ADNI2_70-75_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};
    pathesMCI = {'ADNI2_65-75_F_MCI_nii', 'ADNI2_65-75_M_MCI_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');
    [mciSignals] = connData2signalsFile(base, pathesMCI, 'mci');

    % calculate connectivity
    [cnFCs, meanCNFC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc');
    [adFCs, meanADFC] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc');
    [mciFCs, meanMCIFC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'fc');

    [cnGCs, meanCNGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'gc');
    [adGCs, meanADGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'gc');
    [mciGCs, meanMCIGC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'gc');

    [cnTEs, meanCNTE] = calculateConnectivity(cnSignals, roiNames, 'cn', 'te');
    [adTEs, meanADTE] = calculateConnectivity(adSignals, roiNames, 'ad', 'te');
    [mciTEs, meanMCITE] = calculateConnectivity(mciSignals, roiNames, 'mci', 'te');
end


function [signals, roiNames] = connData2signalsFile(base, pathes, type)
    % constant value
    ROINUM = 132;
    START = 4;
    ReLOCATE = [1	3	5	7	9	11	13	15	17	19	21	23	25	27	29	31	33	35	37	39	41	43	45	47	50	53	58	60	62	64	66	68	70	72	74	76	78	80	82	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	49	52	55	56	57	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	51	54	59	61	63	65	67	69	71	73	75	77	79	81	83	85	87	89	91	93	95	97	99	101	103	105	107	109	111	113	115	117	119	121	123	125	126	127	128	129	130	131	132];


    % init values
    signals = {};
    roiNames2 = cell(1,ROINUM);
    idx = 0;
    
    outfName = ['data/ad-signal-' type '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    % load CONN data files
    for k=1:length(pathes)
        % load experimental signals
        errorCount = 0;
        subjectNum = 1;
        while errorCount < 2
            sbjfile = sprintf('ROI_Subject%03d_Session001.mat', subjectNum);
            expfile = [base pathes{k} '/conn_project01/data/' sbjfile];

            % check file existance
            if ~exist(expfile, 'file')
                errorCount = errorCount + 1;
                continue;
            end
            % load conn ROI signal file
            load(expfile);
            subjectNum = subjectNum + 1;

            seqLen = size(data{1,START},1);
            si = zeros(ROINUM,seqLen);

            for i=1:ROINUM
                si(i,:) = data{1,START+(i-1)}.';
                roiNames2{1,i} = names{1,START+(i-1)};
            end
            % roi relocation
            si2 = si(ReLOCATE,:);
            signals{end+1} = si2;
            roiNames = roiNames2(1, ReLOCATE);
            
            idx = idx + 1;

            % show signals
            figure;
            plot(si.');
            title([type ':' num2str(ROINUM) ' ROI signals (' num2str(idx) ')']);
        end
    end
    save(outfName, 'signals', 'roiNames');
end

function [weights, meanWeight] = calculateConnectivity(signals, roiNames, type, algorithm)
    % constant value
    ROINUM = 132;
    LAG = 3;
    
    outfName = ['data/ad-' algorithm '-' type '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        weights = zeros(ROINUM, ROINUM, length(signals));
        for i=1:length(signals)
            switch(algorithm)
            case 'fc'
                mat = calcFunctionalConnectivity(signals{i});
            case 'gc'
                mat = calcMultivariateGCI(signals{i}, LAG);
            case 'te'
                mat = calcLinueTE(signals{i}, LAG);
            end
            weights(:,:,i) = mat;
        end
        save(outfName, 'weights', 'roiNames');
    end
    meanWeight = nanmean(weights, 3);

    % show functional conectivity
    figure; 
    switch(algorithm)
    case 'fc'
        clims = [-1,1];
        titleStr = [type ' : Functional Connectivity'];
    case 'gc'
        sigma = std(meanWeight(:),'omitnan');
        avg = mean(meanWeight(:),'omitnan');
        meanWeight = (meanWeight - avg) / sigma;
        clims = [-3, 3];
        titleStr = [type ' : multivariate Granger Causality Index'];
    case 'te'
        sigma = std(meanWeight(:),'omitnan');
        avg = mean(meanWeight(:),'omitnan');
        meanWeight = (meanWeight - avg) / sigma;
        clims = [-5, 5];
        titleStr = [type ' : Transfer Entropy (LINER)'];
    end
    imagesc(meanWeight,clims);
    daspect([1 1 1]);
    title(titleStr);
    colorbar;
end