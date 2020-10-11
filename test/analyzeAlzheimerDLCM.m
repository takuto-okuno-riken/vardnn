% Before using this function, download Dlingam-1.2 codes from
% https://sites.google.com/site/sshimizu06/Dlingamcode
% and add a path "Dlingam-1.2" and sub folders. And also download kernel-ICA 1.2 code from
% https://www.di.ens.fr/~fbach/kernel-ica/index.htm
% and add a path "kernel-ica1_2" and sub folders.

function analyzeAlzheimerDLCM
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};
    pathesMCI = {'ADNI2_65-75_F_MCI_nii', 'ADNI2_65-75_M_MCI_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');
    [mciSignals] = connData2signalsFile(base, pathesMCI, 'mci');

    % calculate connectivity
    algNum = 12;
    [cnFCs, meanCNFC, stdCNFC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc');
    [adFCs, meanADFC, stdADFC] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc');
    [mciFCs, meanMCIFC, stdMCIFC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'fc');

    [cnPCs, meanCNPC, stdCNPC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pc');
    [adPCs, meanADPC, stdADPC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pc');
    [mciPCs, meanMCIPC, stdMCIPC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'pc');

    [cnWCSs, meanCNWCS, stdCNWCS] = calculateConnectivity(cnSignals, roiNames, 'cn', 'wcs');
    [adWCSs, meanADWCS, stdADWCS] = calculateConnectivity(adSignals, roiNames, 'ad', 'wcs');
    [mciWCSs, meanMCIWCS, stdMCIWCS] = calculateConnectivity(mciSignals, roiNames, 'mci', 'wcs');

    [cnGCs, meanCNGC, stdCNGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'gc');
    [adGCs, meanADGC, stdADGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'gc');
    [mciGCs, meanMCIGC, stdMCIGC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'gc');

    [cnPGCs, meanCNPGC, stdCNPGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pgc');
    [adPGCs, meanADPGC, stdADPGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pgc');
    [mciPGCs, meanMCIPGC, stdMCIPGC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'pgc');

    [cnTEs, meanCNTE, stdCNTE] = calculateConnectivity(cnSignals, roiNames, 'cn', 'te');
    [adTEs, meanADTE, stdADTE] = calculateConnectivity(adSignals, roiNames, 'ad', 'te');
    [mciTEs, meanMCITE, stdMCITE] = calculateConnectivity(mciSignals, roiNames, 'mci', 'te');

    [cnDLs, meanCNDL, stdCNDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm');
    [adDLs, meanADDL, stdADDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm');
    [mciDLs, meanMCIDL, stdMCIDL] = calculateConnectivity(mciSignals, roiNames, 'mci', 'dlcm');

    [cnDLWs, meanCNDLW, stdCNDLW] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlw');
    [adDLWs, meanADDLW, stdADDLW] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlw');
    [mciDLWs, meanMCIDLW, stdMCIDLW] = calculateConnectivity(mciSignals, roiNames, 'mci', 'dlw');

    [cnDLGs, meanCNDLG, stdCNDLG] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlg');
    [adDLGs, meanADDLG, stdADDLG] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlg');
    [mciDLGs, meanMCIDLG, stdMCIDLG] = calculateConnectivity(mciSignals, roiNames, 'mci', 'dlg');

    [cnPCSs, meanCNPCS, stdCNPCS] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pcs');
    [adPCSs, meanADPCS, stdADPCS] = calculateConnectivity(adSignals, roiNames, 'ad', 'pcs');
    [mciPCSs, meanMCIPCS, stdMCIPCS] = calculateConnectivity(mciSignals, roiNames, 'mci', 'pcs');

    [cnCPCs, meanCNCPC, stdCNCPC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'cpc');
    [adCPCs, meanADCPC, stdADCPC] = calculateConnectivity(adSignals, roiNames, 'ad', 'cpc');
    [mciCPCs, meanMCICPC, stdMCICPC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'cpc');

    [cnFGESs, meanCNFGES, stdCNFGES] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fges');
    [adFGESs, meanADFGES, stdADFGES] = calculateConnectivity(adSignals, roiNames, 'ad', 'fges');
    [mciFGESs, meanMCIFGES, stdMCIFGES] = calculateConnectivity(mciSignals, roiNames, 'mci', 'fges');

    % plot correlation and cos similarity
    nanx = eye(size(meanADFC,1),size(meanADFC,2));
    nanx(nanx==1) = NaN;
%    figure; cnadFCr = plotTwoSignalsCorrelation(meanCNFC+nanx, meanADFC+nanx);
%    figure; cnadGCr = plotTwoSignalsCorrelation(meanCNGC, meanADGC);
%    figure; cnadTEr = plotTwoSignalsCorrelation(meanCNTE, meanADTE);
%    figure; cnadDLr = plotTwoSignalsCorrelation(meanCNDL, meanADDL);
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCNFC+nanx, meanADFC+nanx);
    cosSim(2) = getCosSimilarity(meanCNPC+nanx, meanADPC+nanx);
    cosSim(3) = getCosSimilarity(meanCNWCS+nanx, meanADWCS+nanx);
    cosSim(4) = getCosSimilarity(meanCNGC, meanADGC);
    cosSim(5) = getCosSimilarity(meanCNPGC, meanADPGC);
    cosSim(6) = getCosSimilarity(meanCNTE, meanADTE);
    cosSim(7) = getCosSimilarity(meanCNDL, meanADDL);
    cosSim(8) = getCosSimilarity(meanCNDLW, meanADDLW);
    cosSim(9) = getCosSimilarity(meanCNDLG, meanADDLG);
    cosSim(10) = getCosSimilarity(meanCNPCS+nanx, meanADPCS+nanx);
    cosSim(11) = getCosSimilarity(meanCNCPC+nanx, meanADCPC+nanx);
    cosSim(12) = getCosSimilarity(meanCNFGES+nanx, meanADFGES+nanx);
    X = categorical({'FC','PC','WCS','GC','PGC','TE','DLCM-GC','DLW','dLiNG','PCS','CPC','FGES'});
    figure; bar(X, cosSim);
    title('cos similarity between CN and AD by each algorithm');
    
    % normality test
%{
    cnFCsNt = calculateAlzNormalityTest(cnFCs, roiNames, 'cn', 'fc');
    adFCsNt = calculateAlzNormalityTest(adFCs, roiNames, 'ad', 'fc');
    mciFCsNt = calculateAlzNormalityTest(mciFCs, roiNames, 'mci', 'fc');

    cnPCsNt = calculateAlzNormalityTest(cnPCs, roiNames, 'cn', 'pc');
    adPCsNt = calculateAlzNormalityTest(adPCs, roiNames, 'ad', 'pc');
    mciPCsNt = calculateAlzNormalityTest(mciPCs, roiNames, 'mci', 'pc');

    cnWCSsNt = calculateAlzNormalityTest(cnWCSs, roiNames, 'cn', 'wcs');
    adWCSsNt = calculateAlzNormalityTest(adWCSs, roiNames, 'ad', 'wcs');
    mciWCSsNt = calculateAlzNormalityTest(mciWCSs, roiNames, 'mci', 'wcs');

    cnGCsNt = calculateAlzNormalityTest(cnGCs, roiNames, 'cn', 'gc');
    adGCsNt = calculateAlzNormalityTest(adGCs, roiNames, 'ad', 'gc');
    mciGCsNt = calculateAlzNormalityTest(mciGCs, roiNames, 'mci', 'gc');

    cnPGCsNt = calculateAlzNormalityTest(cnPGCs, roiNames, 'cn', 'pgc');
    adPGCsNt = calculateAlzNormalityTest(adPGCs, roiNames, 'ad', 'pgc');
    mciPGCsNt = calculateAlzNormalityTest(mciPGCs, roiNames, 'mci', 'pgc');

    cnTEsNt = calculateAlzNormalityTest(cnTEs, roiNames, 'cn', 'te');
    adTEsNt = calculateAlzNormalityTest(adTEs, roiNames, 'ad', 'te');
    mciTEsNt = calculateAlzNormalityTest(mciTEs, roiNames, 'mci', 'te');

    cnDLsNt = calculateAlzNormalityTest(cnDLs, roiNames, 'cn', 'dlcm');
    adDLsNt = calculateAlzNormalityTest(adDLs, roiNames, 'ad', 'dlcm');
    mciDLsNt = calculateAlzNormalityTest(mciDLs, roiNames, 'mci', 'dlcm');

    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cn', 'dlw');
    adDLWsNt = calculateAlzNormalityTest(adDLWs, roiNames, 'ad', 'dlw');
    mciDLWsNt = calculateAlzNormalityTest(mciDLWs, roiNames, 'mci', 'dlw');

    cnDLGsNt = calculateAlzNormalityTest(cnDLGs, roiNames, 'cn', 'dlg');
    adDLGsNt = calculateAlzNormalityTest(adDLGs, roiNames, 'ad', 'dlg');
    mciDLGsNt = calculateAlzNormalityTest(mciDLGs, roiNames, 'mci', 'dlg');

    cnPCSsNt = calculateAlzNormalityTest(cnPCSs, roiNames, 'cn', 'pcs');
    adPCSsNt = calculateAlzNormalityTest(adPCSs, roiNames, 'ad', 'pcs');
    mciPCSsNt = calculateAlzNormalityTest(mciPCSs, roiNames, 'mci', 'pcs');

    cnCPCsNt = calculateAlzNormalityTest(cnCPCs, roiNames, 'cn', 'cpc');
    adCPCsNt = calculateAlzNormalityTest(adCPCs, roiNames, 'ad', 'cpc');
    mciCPCsNt = calculateAlzNormalityTest(mciCPCs, roiNames, 'mci', 'cpc');

    cnFGESsNt = calculateAlzNormalityTest(cnFGESs, roiNames, 'cn', 'fges');
    adFGESsNt = calculateAlzNormalityTest(adFGESs, roiNames, 'ad', 'fges');
    mciFGESsNt = calculateAlzNormalityTest(mciFGESs, roiNames, 'mci', 'fges');
%}
    % compalizon test (Wilcoxon, Mann?Whitney U test)
    [cnadFCsUt, cnadFCsUtP, cnadFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, adFCs, roiNames, 'cn', 'ad', 'fc');
    [cnadPCsUt, cnadPCsUtP, cnadPCsUtP2] = calculateAlzWilcoxonTest(cnPCs, adPCs, roiNames, 'cn', 'ad', 'pc');
    [cnadWCSsUt, cnadWCSsUtP, cnadWCSsUtP2] = calculateAlzWilcoxonTest(cnWCSs, adWCSs, roiNames, 'cn', 'ad', 'wcs');
    [cnadGCsUt, cnadGCsUtP, cnadGCsUtP2] = calculateAlzWilcoxonTest(cnGCs, adGCs, roiNames, 'cn', 'ad', 'gc');
    [cnadPGCsUt, cnadPGCsUtP, cnadPGCsUtP2] = calculateAlzWilcoxonTest(cnPGCs, adPGCs, roiNames, 'cn', 'ad', 'pgc');
    [cnadTEsUt, cnadTEsUtP, cnadTEsUtP2] = calculateAlzWilcoxonTest(cnTEs, adTEs, roiNames, 'cn', 'ad', 'te');
    [cnadDLsUt, cnadDLsUtP, cnadDLsUtP2] = calculateAlzWilcoxonTest(cnDLs, adDLs, roiNames, 'cn', 'ad', 'dlcm');
    [cnadDLWsUt, cnadDLWsUtP, cnadDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, adDLWs, roiNames, 'cn', 'ad', 'dlw');
    [cnadDLGsUt, cnadDLGsUtP, cnadDLGsUtP2] = calculateAlzWilcoxonTest(cnDLGs, adDLGs, roiNames, 'cn', 'ad', 'dlg');
    [cnadPCSsUt, cnadPCSsUtP, cnadPCSsUtP2] = calculateAlzWilcoxonTest(cnPCSs, adPCSs, roiNames, 'cn', 'ad', 'pcs');
    [cnadCPCsUt, cnadCPCsUtP, cnadCPCsUtP2] = calculateAlzWilcoxonTest(cnCPCs, adCPCs, roiNames, 'cn', 'ad', 'cpc');
    [cnadFGESsUt, cnadFGESsUtP, cnadFGESsUtP2] = calculateAlzWilcoxonTest(cnFGESs, adFGESs, roiNames, 'cn', 'ad', 'fges');
    
    % using minimum 100 p-value relations. perform 5-fold cross validation.
    topNum = 100;
    sigTh = 2;
    N = 5;

    fcAUC = zeros(1,N);
    pcAUC = zeros(1,N);
    wcsAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    dlgAUC = zeros(1,N);
    teAUC = zeros(1,N);
    pcsAUC = zeros(1,N);
    cpcAUC = zeros(1,N);
    fgesAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);
    dlgROC = cell(N,2);
    teROC = cell(N,2);
    pcsROC = cell(N,2);
    cpcROC = cell(N,2);
    fgesROC = cell(N,2);

    sigCntCN = cell(N,algNum);
    sigCntAD = cell(N,algNum);
    for k=1:N
        i = 1;
        % check sigma of healthy subject
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnFCs, adFCs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, cnadFCsUtP, topNum);
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [fcROC{k,1}, fcROC{k,2}, fcAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCs, adPCs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, cnadPCsUtP, topNum);
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [pcROC{k,1}, pcROC{k,2}, pcAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnWCSs, adWCSs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, cnadWCSsUtP, topNum);
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [wcsROC{k,1}, wcsROC{k,2}, wcsAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnGCs, adGCs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, cnadGCsUtP, topNum);
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [gcROC{k,1}, gcROC{k,2}, gcAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPGCs, adPGCs, k, N);        % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadPGCsUtP, topNum);                                 % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);        % replace *ROC, *AUC
        
        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnTEs, adTEs, k, N);        % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadTEsUtP, topNum);                                % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [teROC{k,1}, teROC{k,2}, teAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLs, adDLs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, cnadDLsUtP, topNum);
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlROC{k,1}, dlROC{k,2}, dlAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs, adDLWs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadDLWsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLGs, adDLGs, k, N);
        [B, I, X] = sortAndPairPValues(control, target, cnadDLGsUtP, topNum);
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlgROC{k,1}, dlgROC{k,2}, dlgAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);
        
        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCSs, adPCSs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadPCSsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [pcsROC{k,1}, pcsROC{k,2}, pcsAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
%{        
        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnCPCs, adCPCs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadCPCsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [cpcROC{k,1}, cpcROC{k,2}, cpcAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnFGESs, adFGESs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadFGESsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [fgesROC{k,1}, fgesROC{k,2}, fgesAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
%}
    end
    figure; boxplot(X);

    % save result
    fname = ['results/ad-cn-ad-roi' num2str(132) '-result.mat'];
    save(fname, 'cosSim', 'fcAUC','pcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','dlgAUC','teAUC','pcsAUC','cpcAUC','fgesAUC', 'fcROC','pcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','dlgROC','teROC','pcsROC','cpcROC','fgesROC', 'sigCntCN', 'sigCntAD');
    mean(dlAUC) % show result AUC
    mean(dlwAUC) % show result AUC
    mean(fcAUC) % show result AUC
    mean(pgcAUC) % show result AUC
    mean(pcAUC) % show result AUC
    mean(wcsAUC) % show result AUC
    mean(dlgAUC) % show result AUC
    mean(pcsAUC) % show result AUC

    % plot ROC curve
    figure;
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(wcsROC, N, [0.9,0.5,0]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
    plotErrorROCcurve(pgcROC, N, [0.0,0.5,0.0]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % TODO:
%    plotErrorROCcurve(dcmROC, N, [0.2,0.2,0.8]);
%    plotErrorROCcurve(rnnROC, N, [0.8,0.8,0.2]);
    plotErrorROCcurve(teROC, N, [0.2,0.6,0.8]);
%    plotErrorROCcurve(nnnueROC, N, [0.8,0.2,0.8]);
    plotErrorROCcurve(dlgROC, N, [0.6,0.6,0.3]);
    plotErrorROCcurve(pcsROC, N, [0.5,0.5,0.5]);
    plotErrorROCcurve(cpcROC, N, [0.5,0.5,0.5]);
    plotErrorROCcurve(fgesROC, N, [0.5,0.5,0.5]);
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '--', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
    plotAverageROCcurve(pgcROC, N, '--', [0.0,0.5,0.0],0.5);
    plotAverageROCcurve(dlROC, N, '-', [0.2,0.2,0.2],1.2);
    plotAverageROCcurve(dlwROC, N, '--', [0.2,0.2,0.2],0.7); % TODO:
%    plotAverageROCcurve(dcmROC, N, '-', [0.2,0.2,0.8],0.5);
%    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5);
    plotAverageROCcurve(teROC, N, '--', [0.2,0.5,0.7],0.5);
%    plotAverageROCcurve(nnnueROC, N, '--', [0.7,0.2,0.7],0.5);
    plotAverageROCcurve(dlgROC, N, '-.', [0.6,0.6,0.3],0.5);
    plotAverageROCcurve(pcsROC, N, '-', [0.5,0.5,0.5],0.5);
    plotAverageROCcurve(cpcROC, N, '--', [0.5,0.5,0.5],0.5);
    plotAverageROCcurve(fgesROC, N, '-.', [0.5,0.5,0.5],0.5);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve']);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

function [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(orgControl, orgTarget, k, N)
    un = floor(size(orgControl,3) / N);
    st = (k-1)*un+1;
    ed = k*un;
    if k==N, ed = size(orgControl,3); end
    control = orgControl(:,:,st:ed);
    un = floor(size(orgTarget,3) / N);
    st = (k-1)*un+1;
    ed = k*un;
    if k==N, ed = size(orgTarget,3); end
    target = orgTarget(:,:,st:ed);
    if N > 1
        orgTarget(:,:,st:ed) = [];
    end
    meanTarget = nanmean(orgTarget, 3);
    stdTarget = nanstd(orgTarget, 1, 3);
    meanControl = nanmean(orgControl, 3);
end

function [x, y, auc] = invertROCcurve(inx, iny)
    y = inx;
    x = iny;
    auc = trapz(x, y);
end

function [x, y, auc] = calcAlzROCcurve(control, target, start)
    x = [0]; y = [0]; % start from (0,0)
    tpmax = length(control);
    fpmax = length(target);
    for i=start:-1:0
        tp = length(find(control>=i));
        fp = length(find(target>=i));
        x = [x fp/fpmax];
        y = [y tp/tpmax];
    end
    auc = trapz(x, y);
end

function [B, I, X] = sortAndPairPValues(control, target, utestP2, topNum)
    ROINUM = size(control,1);
    [B, I] = sort(utestP2(:));
    X = [];
    for k=1:topNum
        i = floor(mod(I(k)-1,ROINUM) + 1);
        j = floor((I(k)-1)/ROINUM + 1);
        x = squeeze(control(i,j,:));
        y = squeeze(target(i,j,:));
        
        if length(x) > length(y)
            x2 = x;
            y2 = nan(length(x),1);
            y2(1:length(y),1) = y;
        else
            x2 = nan(length(y),1);
            x2(1:length(x),1) = x;
            y2 = y;
        end
        X = [X, x2, y2];
    end
end

function sigCount = calcAlzSigmaSubjects(weights, meanWeights, stdWeights, meanControl, I, topNum, sigTh)
    ROINUM = size(weights,1);
    subjectNum = size(weights,3);
    X = [];
    sigCount = [];
    isControlBig = meanControl - meanWeights;
    for n=1:subjectNum
        w2 = weights(:,:,n);
        w2 = w2 - meanWeights;
        sig = abs(w2 ./ stdWeights);
%        sig = w2 ./ stdWeights;
        s = nan(topNum, 1);
        for k=1:topNum
            i = floor(mod(I(k)-1,ROINUM) + 1);
            j = floor((I(k)-1)/ROINUM + 1);
            s(k) = sig(i, j);
%{
            s2 = sig(i, j);
            if s2 < 0 && isControlBig(i, j) > 0
                s2 = 0;
            elseif s2 > 0 && isControlBig(i, j) < 0
                s2 = 0;
            end
            s(k) = abs(s2);
%}
        end
        X = [X, s];
        sigCount = [sigCount, length(find(s>=sigTh))];
    end
%    figure; boxplot(X);
%    figure; bar(sigCount);
end

function [signals, roiNames] = connData2signalsFile(base, pathes, group)
    % constant value
    ROINUM = 132;
    START = 4;
    ReLOCATE = [1	3	5	7	9	11	13	15	17	19	21	23	25	27	29	31	33	35	37	39	41	43	45	47	50	53	58	60	62	64	66	68	70	72	74	76	78	80	82	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	49	52	55	56	57	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	51	54	59	61	63	65	67	69	71	73	75	77	79	81	83	85	87	89	91	93	95	97	99	101	103	105	107	109	111	113	115	117	119	121	123	125	126	127	128	129	130	131	132];

    % init values
    signals = {};
    roiNames2 = cell(1,ROINUM);
    idx = 0;
    
    outfName = ['data/ad-signal-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    % load CONN data files
    for k=1:length(pathes)
        % load experimental signals
        errorCount = 0;
        subjectNum = 0;
        while errorCount < 3
            subjectNum = subjectNum + 1;
            sbjfile = sprintf('ROI_Subject%03d_Session001.mat', subjectNum);
            expfile = [base pathes{k} '/conn_project01/data/' sbjfile];

            % check file existance
            if ~exist(expfile, 'file')
                errorCount = errorCount + 1;
                continue;
            end
            % load conn ROI signal file
            load(expfile);
            errorCount = 0;

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
            title([group ':' num2str(ROINUM) ' ROI signals (' num2str(idx) ')']);
        end
    end
    save(outfName, 'signals', 'roiNames');
end

function [weights, meanWeights, stdWeights] = calculateConnectivity(signals, roiNames, group, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    LAG = 3;

    outfName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        weights = zeros(ROINUM, ROINUM, length(signals));
        for i=1:length(signals)
            switch(algorithm)
            case 'fc'
                mat = calcFunctionalConnectivity(signals{i});
            case 'pc'
                mat = calcPartialCorrelation(signals{i});
            case 'wcs'
                fName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(fName, 'file')
                    load(fName);
                else
                    mat = calcWaveletCoherence(signals{i});
                    save(fName, 'mat');
                end
            case 'gc'
                mat = calcMultivariateGCI(signals{i}, LAG);
            case 'pgc'
                mat = calcPairwiseGCI(signals{i}, LAG);
            case 'te'
                mat = calcLinueTE(signals{i}, LAG);
            case 'pcs'
                csvFile = ['results/tetrad/pcs-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'cpc'
                csvFile = ['results/tetrad/cpc-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'fges'
                csvFile = ['results/tetrad/fges-ad-signal-' group '-roi' num2str(ROINUM) '-' num2str(i) '.csv'];
                mat = readmatrix(csvFile);
            case 'dlg'
                fName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(fName, 'file')
                    load(fName);
                else
                    mat = calcDirectLiNGAM(signals{i});
                    save(fName, 'mat');
                end
            case 'dlcm'
                dlcmName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                if exist(dlcmName, 'file')
                    load(dlcmName);
                else
                    [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
                    sigLen = size(si,2);
                    inSignal = rand(ROINUM, sigLen);
                    inControl = eye(ROINUM);
                    netDLCM = initDlcmNetwork(si, inSignal, [], inControl); 
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
                    netDLCM = trainDlcmNetwork(si, inSignal, [], inControl, netDLCM, options);
                    [time, loss, rsme] = getDlcmTrainingResult(netDLCM);
                    disp(['end training : rsme=' num2str(rsme)]);
                    % calc dlcm-gc
                    mat = calcDlcmGCI(si, inSignal, [], inControl, netDLCM);
                    
                    save(dlcmName, 'netDLCM', 'si', 'inSignal', 'inControl', 'mat', 'sig', 'c', 'maxsi', 'minsi');
                end
            case 'dlw' % should be called after dlcm
                dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                load(dlcmName);
                mat = calcDlcmWCI(netDLCM, [], inControl);
            end
            weights(:,:,i) = mat;
        end
        save(outfName, 'weights', 'roiNames');
    end
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);
    % counting by source region and target region
    meanSource = nanmean(meanWeights,1);
    meanTarget = nanmean(meanWeights,2);
    save(outfName, 'weights', 'meanWeights', 'stdWeights', 'roiNames', 'meanSource', 'meanTarget');

    % show functional conectivity
    figure; 
    switch(algorithm)
    case 'fc'
        clims = [-1,1];
        titleStr = [group ' : Functional Connectivity'];
        sigWeights = meanWeights;
    case 'pc'
        clims = [-1,1];
        titleStr = [group ' : Partial Correlation'];
        sigWeights = meanWeights;
    case 'wcs'
        clims = [-1,1];
        titleStr = [group ' : Wavelet Coherence Cross Spectrum'];
        sigWeights = meanWeights;
    case 'gc'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : multivariate Granger Causality Index'];
    case 'pgc'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : pairwise Granger Causality Index'];
    case 'te'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : Transfer Entropy (LINER)'];
    case 'pcs'
        clims = [-1,1];
        titleStr = [group ' : PC-stable-max'];
        sigWeights = meanWeights;
    case 'cpc'
        clims = [-1,1];
        titleStr = [group ' : Conservative PC'];
        sigWeights = meanWeights;
    case 'fges'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : FGES'];
    case 'dlg'
        clims = [-1,1];
        titleStr = [group ' : Direct LiNGAM'];
        sigWeights = meanWeights;
    case 'dlcm'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM Granger Causality Index'];
    case 'dlw'
        sigma = std(meanWeights(:),'omitnan');
        avg = mean(meanWeights(:),'omitnan');
        sigWeights = (meanWeights - avg) / sigma;
        clims = [-3, 3];
        titleStr = [group ' : DLCM Weight Causality Index'];
    end
    imagesc(sigWeights,clims);
    daspect([1 1 1]);
    title(titleStr);
    colorbar;
end

function [normalities, normalitiesP] = calculateAlzNormalityTest(weights, roiNames, group, algorithm)
    % constant value
    ROINUM = size(weights,1);

    outfName = ['results/ad-' algorithm '-' group '-roi' num2str(ROINUM) '-normality.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        normalities = nan(ROINUM, ROINUM);
        normalitiesP = nan(ROINUM, ROINUM);
        for i=1:ROINUM
            for j=1:ROINUM
                if i==j, continue; end
                x = squeeze(weights(i,j,:));
                [h, p] = lillietest(x);
                normalities(i,j) = 1 - h;
                normalitiesP(i,j) = p;
            end
        end
        save(outfName, 'normalities', 'normalitiesP', 'roiNames');
    end

    load('test/colormap.mat')
    % show normality test result
    figure; 
    colormap(hvalmap);
    clims = [0,1];
    imagesc(normalities,clims);
    daspect([1 1 1]);
    title([group '-' algorithm ' normality test result']);
    colorbar;
    % normality test p values
    normalitiesP(isnan(normalitiesP)) = 1;
    figure;
    colormap(pvalmap);
    clims = [0,0.5];
    imagesc(normalitiesP,clims);
    daspect([1 1 1]);
    title([group '-' algorithm ' normality test p values']);
    colorbar;
end


function [utestH, utestP, utestP2] = calculateAlzWilcoxonTest(control, target, roiNames, controlGroup, targetGroup, algorithm)
    % constant value
    ROINUM = size(control,1);

    outfName = ['results/ad-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROINUM) '-utest.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        utestH = nan(ROINUM, ROINUM);
        utestP = nan(ROINUM, ROINUM);
        utestP2 = nan(ROINUM, ROINUM);
        for i=1:ROINUM
            for j=1:ROINUM
                if i==j, continue; end
                x = squeeze(control(i,j,:));
                y = squeeze(target(i,j,:));
                if length(x) > length(y)
                    x2 = x;
                    y2 = nan(length(x),1);
                    y2(1:length(y),1) = y;
                else
                    x2 = nan(length(y),1);
                    x2(1:length(x),1) = x;
                    y2 = y;
                end
                if isempty(find(~isnan(x2))) || isempty(find(~isnan(y2))), continue; end
                [p, h] = ranksum(x2,y2);
%                [p, h] = signrank(x2,y2);
                utestH(i,j) = h;
                utestP(i,j) = p;
                if h > 0 && nanmean(x) > nanmean(y)
                    utestP2(i,j) = p;
                end
            end
        end
        save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames');
    end
    % counting by source region and target region
    countSource = nansum(utestH,1);
    countTarget = nansum(utestH,2);
    save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames', 'countSource', 'countTarget');

    load('test/colormap.mat')
    % U test result
    figure; 
    colormap(hvalmap);
    clims = [0,1];
    imagesc(utestH,clims);
    daspect([1 1 1]);
    title([controlGroup '-' targetGroup ' : ' algorithm ' : u test result']);
    colorbar;
    % U test p values
    utestP(isnan(utestP)) = 1;
    figure;
    colormap(pvalmap);
    clims = [0,1];
    imagesc(utestP, clims);
    daspect([1 1 1]);
    title([controlGroup '-' targetGroup ' : ' algorithm ' : u test p values']);
    colorbar;
end
