% Before using this function, download Dlingam-1.2 codes from
% https://sites.google.com/site/sshimizu06/Dlingamcode
% and add a path "Dlingam-1.2" and sub folders. And also download kernel-ICA 1.2 code from
% https://www.di.ens.fr/~fbach/kernel-ica/index.htm
% and add a path "kernel-ica1_2" and sub folders.

% Before using this function, download PartiallyConditionedGrangerCausality codes from
% https://github.com/danielemarinazzo/PartiallyConditionedGrangerCausality
% and add a path "PartiallyConditionedGrangerCausality-master" and sub folders. 

% original alzheimer diagnosis. 5-fold cross validation with top 100 regional relation
% based on across all subject. (may be biased case)

function analyzeAlzheimerVarDnn
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};
    pathesMCI = {'ADNI2_65-75_F_MCI_nii', 'ADNI2_65-75_M_MCI_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn', 'data/ad', 'ad');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad', 'data/ad', 'ad');
    [mciSignals] = connData2signalsFile(base, pathesMCI, 'mci', 'data/ad', 'ad');

    global resultsPath;
    global resultsPrefix;
    resultsPath = 'results/ad';
    resultsPrefix = 'ad';

    % check amplitude difference by ROI
    [ampDiffROI,cnStds,adStds] = checkAmplitudeDiffROI(cnSignals, adSignals);

    % calculate connectivity
    algNum = 31;
    [cnFCs, meanCNFC, stdCNFC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc');
    [adFCs, meanADFC, stdADFC] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc');
    [mciFCs, meanMCIFC, stdMCIFC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'fc');

    [cnPCs, meanCNPC, stdCNPC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pc');
    [adPCs, meanADPC, stdADPC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pc');
    [mciPCs, meanMCIPC, stdMCIPC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'pc');

    [cnPcPCs, meanCNPcPC, stdCNPcPC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pcpc');
    [adPcPCs, meanADPcPC, stdADPcPC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pcpc');

    [cnLsoPCs, meanCNLsoPC, stdCNLsoPC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'lsopc');
    [adLsoPCs, meanADLsoPC, stdADLsoPC] = calculateConnectivity(adSignals, roiNames, 'ad', 'lsopc');

    [cnPlsPCs, meanCNPlsPC, stdCNPlsPC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'plspc');
    [adPlsPCs, meanADPlsPC, stdADPlsPC] = calculateConnectivity(adSignals, roiNames, 'ad', 'plspc');
    
    [cnWCSs, meanCNWCS, stdCNWCS] = calculateConnectivity(cnSignals, roiNames, 'cn', 'wcs');
    [adWCSs, meanADWCS, stdADWCS] = calculateConnectivity(adSignals, roiNames, 'ad', 'wcs');
    [mciWCSs, meanMCIWCS, stdMCIWCS] = calculateConnectivity(mciSignals, roiNames, 'mci', 'wcs');

    [cnGCs, meanCNGC, stdCNGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'gc', 0, 3);
    [adGCs, meanADGC, stdADGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'gc', 0, 3);
    [mciGCs, meanMCIGC, stdMCIGC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'gc', 0, 3);

    [cnPGCs, meanCNPGC, stdCNPGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pgc', 0, 3);
    [adPGCs, meanADPGC, stdADPGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pgc', 0, 3);
    [mciPGCs, meanMCIPGC, stdMCIPGC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'pgc', 0, 3);

    [cnTEs, meanCNTE, stdCNTE] = calculateConnectivity(cnSignals, roiNames, 'cn', 'te', 0, 3);
    [adTEs, meanADTE, stdADTE] = calculateConnectivity(adSignals, roiNames, 'ad', 'te', 0, 3);
    [mciTEs, meanMCITE, stdMCITE] = calculateConnectivity(mciSignals, roiNames, 'mci', 'te', 0, 3);

    [cnDLs, meanCNDL, stdCNDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 0, 1, 1);
    [adDLs, meanADDL, stdADDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 0, 1, 1);
    [mciDLs, meanMCIDL, stdMCIDL] = calculateConnectivity(mciSignals, roiNames, 'mci', 'dlcm', 0, 1, 1);

    [cnDLWs, meanCNDLW, stdCNDLW] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlw', 0, 1, 1);
    [adDLWs, meanADDLW, stdADDLW] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlw', 0, 1, 1);
    [mciDLWs, meanMCIDLW, stdMCIDLW] = calculateConnectivity(mciSignals, roiNames, 'mci', 'dlw', 0, 1, 1);

    [cnPCDLs, meanCNPCDL, stdCNPCDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pcdl', 0, 1, 1);
    [adPCDLs, meanADPCDL, stdADPCDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'pcdl', 0, 1, 1);

    [cnPCDLWs, meanCNPCDLW, stdCNPCDLW] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pcdlw', 0, 1, 1);
    [adPCDLWs, meanADPCDLW, stdADPCDLW] = calculateConnectivity(adSignals, roiNames, 'ad', 'pcdlw', 0, 1, 1);
    
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

    [cnFCAs, meanCNFCA, stdCNFCA] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fca');
    [adFCAs, meanADFCA, stdADFCA] = calculateConnectivity(adSignals, roiNames, 'ad', 'fca');
    [mciFCAs, meanMCIFCA, stdMCIFCA] = calculateConnectivity(mciSignals, roiNames, 'mci', 'fca');

    [cnTSFCs, meanCNTSFC, stdCNTSFC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'tsfc', 0, 3);
    [adTSFCs, meanADTSFC, stdADTSFC] = calculateConnectivity(adSignals, roiNames, 'ad', 'tsfc', 0, 3);
    [mciTSFCs, meanMCITSFC, stdMCITSFC] = calculateConnectivity(mciSignals, roiNames, 'mci', 'tsfc', 0, 3);

    [cnTSFCAs, meanCNTSFCA, stdCNTSFCA] = calculateConnectivity(cnSignals, roiNames, 'cn', 'tsfca', 0, 3);
    [adTSFCAs, meanADTSFCA, stdADTSFCA] = calculateConnectivity(adSignals, roiNames, 'ad', 'tsfca', 0, 3);
    [mciTSFCAs, meanMCITSFCA, stdMCITSFCA] = calculateConnectivity(mciSignals, roiNames, 'mci', 'tsfca', 0, 3);

    [cnMVARDIs, meanCNMVARDI, stdCNMVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mvarec', 0, 3);
    [adMVARDIs, meanADMVARDI, stdADMVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'mvarec', 0, 3);
    [mciMVARDIs, meanMCIMVARDI, stdMCIMVARDI] = calculateConnectivity(mciSignals, roiNames, 'mci', 'mvarec', 0, 3);

    [cnPVARDIs, meanCNPVARDI, stdCNPVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pvarec', 0, 3);
    [adPVARDIs, meanADPVARDI, stdADPVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'pvarec', 0, 3);

    [cnMPCVARDIs, meanCNMPCVARDI, stdCNMPCVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mpcvarec', 0, 3);
    [adMPCVARDIs, meanADMPCVARDI, stdADMPCVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'mpcvarec', 0, 3);
    [mciMPCVARDIs, meanMCIMPCVARDI, stdMCIMPCVARDI] = calculateConnectivity(mciSignals, roiNames, 'mci', 'mpcvarec', 0, 3);

    [cnMPCVARGCs, meanCNMPCVARGC, stdCNMPCVARGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mpcvargc', 0, 3);
    [adMPCVARGCs, meanADMPCVARGC, stdADMPCVARGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'mpcvargc', 0, 3);

    [cnPPCVARDIs, meanCNPPCVARDI, stdCNPPCVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'ppcvarec', 0, 3);
    [adPPCVARDIs, meanADPPCVARDI, stdADPPCVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'ppcvarec', 0, 3);

    [cnPPCVARGCs, meanCNPPCVARGC, stdCNPPCVARGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'ppcvargc', 0, 3);
    [adPPCVARGCs, meanADPPCVARGC, stdADPPCVARGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'ppcvargc', 0, 3);

    [cnMPLSVARDIs, meanCNMPLSVARDI, stdCNMPLSVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mplsvarec', 0, 3);
    [adMPLSVARDIs, meanADMPLSVARDI, stdADMPLSVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'mplsvarec', 0, 3);

    [cnMPLSVARGCs, meanCNMPLSVARGC, stdCNMPLSVARGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mplsvargc', 0, 3);
    [adMPLSVARGCs, meanADMPLSVARGC, stdADMPLSVARGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'mplsvargc', 0, 3);

    [cnPPLSVARDIs, meanCNPPLSVARDI, stdCNPPLSVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pplsvarec', 0, 3);
    [adPPLSVARDIs, meanADPPLSVARDI, stdADPPLSVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'pplsvarec', 0, 3);

    [cnPPLSVARGCs, meanCNPPLSVARGC, stdCNPPLSVARGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pplsvargc', 0, 3);
    [adPPLSVARGCs, meanADPPLSVARGC, stdADPPLSVARGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pplsvargc', 0, 3);

    [cnMLSOVARDIs, meanCNMLSOVARDI, stdCNMLSOVARDI] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mlsovarec', 0, 3);
    [adMLSOVARDIs, meanADMLSOVARDI, stdADMLSOVARDI] = calculateConnectivity(adSignals, roiNames, 'ad', 'mlsovarec', 0, 3);

    [cnMLSOVARGCs, meanCNMLSOVARGC, stdCNMLSOVARGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'mlsovargc', 0, 3);
    [adMLSOVARGCs, meanADMLSOVARGC, stdADMLSOVARGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'mlsovargc', 0, 3);

    [cnPCGCs, meanCNPCGC, stdCNPCGC] = calculateConnectivity(cnSignals, roiNames, 'cn', 'pcgc', 0, 3);
    [adPCGCs, meanADPCGC, stdADPCGC] = calculateConnectivity(adSignals, roiNames, 'ad', 'pcgc', 0, 3);

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
    cosSim(3) = getCosSimilarity(meanCNPcPC+nanx, meanADPcPC+nanx);
    cosSim(4) = getCosSimilarity(meanCNLsoPC+nanx, meanADLsoPC+nanx);
    cosSim(5) = getCosSimilarity(meanCNPlsPC+nanx, meanADPlsPC+nanx);
    cosSim(6) = getCosSimilarity(meanCNPGC, meanADPGC);
    cosSim(7) = getCosSimilarity(meanCNGC, meanADGC);
    cosSim(8) = getCosSimilarity(meanCNMPCVARGC+nanx, meanADMPCVARGC+nanx);
    cosSim(9) = getCosSimilarity(meanCNMLSOVARGC+nanx, meanADMLSOVARGC+nanx);
    cosSim(10) = getCosSimilarity(meanCNMPLSVARGC+nanx, meanADMPLSVARGC+nanx);
    cosSim(11) = 0; % RNN-GC
    cosSim(12) = getCosSimilarity(meanCNPCGC+nanx, meanADPCGC+nanx);
    cosSim(13) = getCosSimilarity(meanCNTE, meanADTE); % LINUE-TE
    cosSim(14) = 0; % NNNUE-TE
    cosSim(15) = getCosSimilarity(meanCNDL, meanADDL); % VARDNN-GC
    cosSim(16) = getCosSimilarity(meanCNDLW, meanADDLW);
    cosSim(17) = getCosSimilarity(meanCNPCDL, meanADPCDL); % VARDNN-GC
    cosSim(18) = getCosSimilarity(meanCNPCDLW, meanADPCDLW); i=19; 
    cosSim(i) = getCosSimilarity(meanCNWCS+nanx, meanADWCS+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNDLG, meanADDLG); i=i+1; % dLINGAM
    cosSim(i) = getCosSimilarity(meanCNPCS+nanx, meanADPCS+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNCPC+nanx, meanADCPC+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNFGES+nanx, meanADFGES+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNFCA+nanx, meanADFCA+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNTSFC+nanx, meanADTSFC+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNTSFCA+nanx, meanADTSFCA+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNMVARDI+nanx, meanADMVARDI+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNPVARDI+nanx, meanADPVARDI+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNMPCVARDI+nanx, meanADMPCVARDI+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNPPCVARDI+nanx, meanADPPCVARDI+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNPPCVARGC+nanx, meanADPPCVARGC+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNMPLSVARDI+nanx, meanADMPLSVARDI+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNPPLSVARDI+nanx, meanADPPLSVARDI+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNPPLSVARGC+nanx, meanADPPLSVARGC+nanx); i=i+1;
    cosSim(i) = getCosSimilarity(meanCNMLSOVARDI+nanx, meanADMLSOVARDI+nanx); i=i+1;
    figure; bar(cosSim);
    title('cos similarity between CN and AD by each algorithm');
    
    % normality test
%{
    cnFCsNt = calculateAlzNormalityTest(cnFCs, roiNames, 'cn', 'fc');
    adFCsNt = calculateAlzNormalityTest(adFCs, roiNames, 'ad', 'fc');
    mciFCsNt = calculateAlzNormalityTest(mciFCs, roiNames, 'mci', 'fc');

%    cnPCsNt = calculateAlzNormalityTest(cnPCs, roiNames, 'cn', 'pc');
%    adPCsNt = calculateAlzNormalityTest(adPCs, roiNames, 'ad', 'pc');
%    mciPCsNt = calculateAlzNormalityTest(mciPCs, roiNames, 'mci', 'pc');

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
%}
    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cn', 'dlw');
    adDLWsNt = calculateAlzNormalityTest(adDLWs, roiNames, 'ad', 'dlw');
    mciDLWsNt = calculateAlzNormalityTest(mciDLWs, roiNames, 'mci', 'dlw');
%{
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
    [cnadPcPCsUt, cnadPcPCsUtP, cnadPcPCsUtP2] = calculateAlzWilcoxonTest(cnPcPCs, adPcPCs, roiNames, 'cn', 'ad', 'pcpc');
    [cnadLsoPCsUt, cnadLsoPCsUtP, cnadLsoPCsUtP2] = calculateAlzWilcoxonTest(cnLsoPCs, adLsoPCs, roiNames, 'cn', 'ad', 'lsopc');
    [cnadPlsPCsUt, cnadPlsPCsUtP, cnadPlsPCsUtP2] = calculateAlzWilcoxonTest(cnPlsPCs, adPlsPCs, roiNames, 'cn', 'ad', 'plspc');
    [cnadWCSsUt, cnadWCSsUtP, cnadWCSsUtP2] = calculateAlzWilcoxonTest(cnWCSs, adWCSs, roiNames, 'cn', 'ad', 'wcs');
    [cnadGCsUt, cnadGCsUtP, cnadGCsUtP2] = calculateAlzWilcoxonTest(cnGCs, adGCs, roiNames, 'cn', 'ad', 'gc');
    [cnadPGCsUt, cnadPGCsUtP, cnadPGCsUtP2] = calculateAlzWilcoxonTest(cnPGCs, adPGCs, roiNames, 'cn', 'ad', 'pgc');
    [cnadTEsUt, cnadTEsUtP, cnadTEsUtP2] = calculateAlzWilcoxonTest(cnTEs, adTEs, roiNames, 'cn', 'ad', 'te');
    [cnadDLsUt, cnadDLsUtP, cnadDLsUtP2] = calculateAlzWilcoxonTest(cnDLs, adDLs, roiNames, 'cn', 'ad', 'dlcm');
    [cnadDLWsUt, cnadDLWsUtP, cnadDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, adDLWs, roiNames, 'cn', 'ad', 'dlw');
    [cnadPCDLsUt, cnadPCDLsUtP, cnadPCDLsUtP2] = calculateAlzWilcoxonTest(cnPCDLs, adPCDLs, roiNames, 'cn', 'ad', 'pcdl');
    [cnadPCDLWsUt, cnadPCDLWsUtP, cnadPCDLWsUtP2] = calculateAlzWilcoxonTest(cnPCDLWs, adPCDLWs, roiNames, 'cn', 'ad', 'pcdlw');
    [cnadDLGsUt, cnadDLGsUtP, cnadDLGsUtP2] = calculateAlzWilcoxonTest(cnDLGs, adDLGs, roiNames, 'cn', 'ad', 'dlg');
    [cnadPCSsUt, cnadPCSsUtP, cnadPCSsUtP2] = calculateAlzWilcoxonTest(cnPCSs, adPCSs, roiNames, 'cn', 'ad', 'pcs');
    [cnadCPCsUt, cnadCPCsUtP, cnadCPCsUtP2] = calculateAlzWilcoxonTest(cnCPCs, adCPCs, roiNames, 'cn', 'ad', 'cpc');
    [cnadFGESsUt, cnadFGESsUtP, cnadFGESsUtP2] = calculateAlzWilcoxonTest(cnFGESs, adFGESs, roiNames, 'cn', 'ad', 'fges');
    [cnadFCAsUt, cnadFCAsUtP, cnadFCAsUtP2] = calculateAlzWilcoxonTest(cnFCAs, adFCAs, roiNames, 'cn', 'ad', 'fca');
    [cnadTsFCsUt, cnadTsFCsUtP, cnadTsFCsUtP2] = calculateAlzWilcoxonTest(cnTSFCs, adTSFCs, roiNames, 'cn', 'ad', 'tsfc');
    [cnadTsFCAsUt, cnadTsFCAsUtP, cnadTsFCAsUtP2] = calculateAlzWilcoxonTest(cnTSFCAs, adTSFCAs, roiNames, 'cn', 'ad', 'tsfca');
    [cnadMvarDIsUt, cnadMvarDIsUtP, cnadMvarDIsUtP2] = calculateAlzWilcoxonTest(cnMVARDIs, adMVARDIs, roiNames, 'cn', 'ad', 'mvarec');
    [cnadPvarDIsUt, cnadPvarDIsUtP, cnadPvarDIsUtP2] = calculateAlzWilcoxonTest(cnPVARDIs, adPVARDIs, roiNames, 'cn', 'ad', 'pvarec');
    [cnadMpcvarDIsUt, cnadMpcvarDIsUtP, cnadMpcvarDIsUtP2] = calculateAlzWilcoxonTest(cnMPCVARDIs, adMPCVARDIs, roiNames, 'cn', 'ad', 'mpcvarec');
    [cnadMpcvarGCsUt, cnadMpcvarGCsUtP, cnadMpcvarGCsUtP2] = calculateAlzWilcoxonTest(cnMPCVARGCs, adMPCVARGCs, roiNames, 'cn', 'ad', 'mpcvargc');
    [cnadPpcvarDIsUt, cnadPpcvarDIsUtP, cnadPpcvarDIsUtP2] = calculateAlzWilcoxonTest(cnPPCVARDIs, adPPCVARDIs, roiNames, 'cn', 'ad', 'ppcvarec');
    [cnadPpcvarGCsUt, cnadPpcvarGCsUtP, cnadPpcvarGCsUtP2] = calculateAlzWilcoxonTest(cnPPCVARGCs, adPPCVARGCs, roiNames, 'cn', 'ad', 'ppcvargc');
    [cnadMplsvarDIsUt, cnadMplsvarDIsUtP, cnadMplsvarDIsUtP2] = calculateAlzWilcoxonTest(cnMPLSVARDIs, adMPLSVARDIs, roiNames, 'cn', 'ad', 'mplsvarec');
    [cnadMplsvarGCsUt, cnadMplsvarGCsUtP, cnadMplsvarGCsUtP2] = calculateAlzWilcoxonTest(cnMPLSVARGCs, adMPLSVARGCs, roiNames, 'cn', 'ad', 'mplsvargc');
    [cnadPplsvarDIsUt, cnadPplsvarDIsUtP, cnadPplsvarDIsUtP2] = calculateAlzWilcoxonTest(cnPPLSVARDIs, adPPLSVARDIs, roiNames, 'cn', 'ad', 'pplsvarec');
    [cnadPplsvarGCsUt, cnadPplsvarGCsUtP, cnadPplsvarGCsUtP2] = calculateAlzWilcoxonTest(cnPPLSVARGCs, adPPLSVARGCs, roiNames, 'cn', 'ad', 'pplsvargc');
    [cnadMlsovarDIsUt, cnadMlsovarDIsUtP, cnadMlsovarDIsUtP2] = calculateAlzWilcoxonTest(cnMLSOVARDIs, adMLSOVARDIs, roiNames, 'cn', 'ad', 'mlsovarec');
    [cnadMlsovarGCsUt, cnadMlsovarGCsUtP, cnadMlsovarGCsUtP2] = calculateAlzWilcoxonTest(cnMLSOVARGCs, adMLSOVARGCs, roiNames, 'cn', 'ad', 'mlsovargc');
    [cnadPCGCsUt, cnadPCGCsUtP, cnadPCGCsUtP2] = calculateAlzWilcoxonTest(cnPCGCs, adPCGCs, roiNames, 'cn', 'ad', 'pcgc');
    
    % show top 100 most different relations
    topNum = 100;
    [B,I]=sort(cnadDLWsUtP(:));
    figure; semilogy(B(1:topNum)); title('VARDNN(1)');
    figure; hold on; plot(meanCNDLW(I(1:topNum))); plot(meanADDLW(I(1:topNum))); hold off; title('VARDNN(1)');
    
    [B,I]=sort(cnadDLsUtP(:));
    figure; semilogy(B(1:topNum)); title('VARDNN(1)-GC');
    figure; hold on; plot(meanCNDL(I(1:topNum))); plot(meanADDL(I(1:topNum))); hold off; title('VARDNN(1)-GC');
    
    [B,I]=sort(cnadFCAsUtP(:));
    figure; semilogy(B(1:topNum)); title('FC');
    figure; hold on; plot(meanCNFC(I(1:topNum))); plot(meanADFC(I(1:topNum))); hold off; title('FC');
    
    [B,I]=sort(cnadMplsvarGCsUtP(:));
    figure; semilogy(B(1:topNum)); title('mPLS-GC(3)');
    figure; hold on; plot(meanCNMPLSVARGC(I(1:topNum))); plot(meanADMPLSVARGC(I(1:topNum))); hold off; title('mPLS-GC(3)');
    
    % show top3 most different DI histgram
    for i=1:3
        [B,I]=sort(cnadDLWsUtP(:));
        [c, r] = index2rowCol(I(i),132); eg=[0:0.0125:0.4];
        figure; hold on; histogram(cnDLWs(r,c,:),eg); histogram(adDLWs(r,c,:),eg); hold off;
        title(['VARDNN(1)-DI top' num2str(i) ' (' num2str(r) ',' num2str(c) ')']);

        [B,I]=sort(cnadDLsUtP(:));
        [c, r] = index2rowCol(I(i),132); eg=[0:0.0125:0.4];
        figure; hold on; histogram(cnDLs(r,c,:),eg); histogram(adDLs(r,c,:),eg); hold off;
        title(['VARDNN(1)-GC top' num2str(i) ' (' num2str(r) ',' num2str(c) ')']);

        [B,I]=sort(cnadFCAsUtP(:));
        [c, r] = index2rowCol(I(i),132); eg=[-0.2:0.05:1];
        figure; hold on; histogram(cnFCs(r,c,:),eg); histogram(adFCs(r,c,:),eg); hold off;
        title(['FC top' num2str(i) ' (' num2str(r) ',' num2str(c) ')']);

        [B,I]=sort(cnadMplsvarGCsUtP(:));
        [c, r] = index2rowCol(I(i),132); eg=[3:0.25:9];
        figure; hold on; histogram(cnMPLSVARGCs(r,c,:),eg); histogram(adMPLSVARGCs(r,c,:),eg); hold off;
        title(['mPLS-GC(3) top' num2str(i) ' (' num2str(r) ',' num2str(c) ')']);
    end
    
    % using minimum [30-300] p-value relations.
    sigThs = [1.5:0.1:2.0];
    rrNums = [30:30:300];
    N = length(sigThs)*length(rrNums);

    fcAUC = zeros(1,N);
    pcAUC = zeros(1,N);
    pcpcAUC = zeros(1,N);
    lsopcAUC = zeros(1,N);
    plspcAUC = zeros(1,N);
    wcsAUC = zeros(1,N);
    gcAUC = zeros(1,N);
    pgcAUC = zeros(1,N);
    dlAUC = zeros(1,N);
    dlwAUC = zeros(1,N);
    pcdlAUC = zeros(1,N);
    pcdlwAUC = zeros(1,N);
    dlgAUC = zeros(1,N);
    teAUC = zeros(1,N);
    pcsAUC = zeros(1,N);
    cpcAUC = zeros(1,N);
    fgesAUC = zeros(1,N);
    fcaAUC = zeros(1,N);
    tsfcAUC = zeros(1,N);
    tsfcaAUC = zeros(1,N);
    mvardiAUC = zeros(1,N);
    pvardiAUC = zeros(1,N);
    mpcvardiAUC = zeros(1,N);
    mpcvargcAUC = zeros(1,N);
    ppcvardiAUC = zeros(1,N);
    ppcvargcAUC = zeros(1,N);
    mplsdiAUC = zeros(1,N);
    mplsgcAUC = zeros(1,N);
    pplsdiAUC = zeros(1,N);
    pplsgcAUC = zeros(1,N);
    mlsodiAUC = zeros(1,N);
    mlsogcAUC = zeros(1,N);
    pcgcAUC = zeros(1,N);
    fcROC = cell(N,2);
    pcROC = cell(N,2);
    pcpcROC = cell(N,2);
    lsopcROC = cell(N,2);
    plspcROC = cell(N,2);
    wcsROC = cell(N,2);
    gcROC = cell(N,2);
    pgcROC = cell(N,2);
    dlROC = cell(N,2);
    dlwROC = cell(N,2);
    pcdlROC = cell(N,2);
    pcdlwROC = cell(N,2);
    dlgROC = cell(N,2);
    teROC = cell(N,2);
    pcsROC = cell(N,2);
    cpcROC = cell(N,2);
    fgesROC = cell(N,2);
    fcaROC = cell(N,2);
    tsfcROC = cell(N,2);
    tsfcaROC = cell(N,2);
    mvardiROC = cell(N,2);
    pvardiROC = cell(N,2);
    mpcvardiROC = cell(N,2);
    mpcvargcROC = cell(N,2);
    ppcvardiROC = cell(N,2);
    ppcvargcROC = cell(N,2);
    mplsdiROC = cell(N,2);
    mplsgcROC = cell(N,2);
    pplsdiROC = cell(N,2);
    pplsgcROC = cell(N,2);
    mlsodiROC = cell(N,2);
    mlsogcROC = cell(N,2);
    pcgcROC = cell(N,2);
    fcACC = cell(N,1);
    pcACC = cell(N,1);
    pcpcACC = cell(N,1);
    lsopcACC = cell(N,1);
    plspcACC = cell(N,1);
    wcsACC = cell(N,1);
    gcACC = cell(N,1);
    pgcACC = cell(N,1);
    dlACC = cell(N,1);
    dlwACC = cell(N,1);
    pcdlACC = cell(N,1);
    pcdlwACC = cell(N,1);
    dlgACC = cell(N,1);
    teACC = cell(N,1);
    pcsACC = cell(N,1);
    cpcACC = cell(N,1);
    fgesACC = cell(N,1);
    fcaACC = cell(N,1);
    tsfcACC = cell(N,1);
    tsfcaACC = cell(N,1);
    mvardiACC = cell(N,1);
    pvardiACC = cell(N,1);
    mpcvardiACC = cell(N,1);
    mpcvargcACC = cell(N,1);
    ppcvardiACC = cell(N,1);
    ppcvargcACC = cell(N,1);
    mplsdiACC = cell(N,1);
    mplsgcACC = cell(N,1);
    pplsdiACC = cell(N,1);
    pplsgcACC = cell(N,1);
    mlsodiACC = cell(N,1);
    mlsogcACC = cell(N,1);
    pcgcACC = cell(N,1);

    nodeNum = size(cnFCs,1);
    pvList = nan(nodeNum*nodeNum,algNum);
    sigCntCN = cell(N,algNum);
    sigCntAD = cell(N,algNum);
    for p=1:length(rrNums)
        topNum = rrNums(p);
        for q=1:length(sigThs)
            sigTh = sigThs(q);

            k = (p-1) * length(sigThs) + q;
            i = 1;
            % check sigma of healthy subject
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnFCs, adFCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadFCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [fcROC{k,1}, fcROC{k,2}, fcAUC(k), fcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCs, adPCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pcROC{k,1}, pcROC{k,2}, pcAUC(k), pcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPcPCs, adPcPCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPcPCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pcpcROC{k,1}, pcpcROC{k,2}, pcpcAUC(k), pcpcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnLsoPCs, adLsoPCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadLsoPCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [lsopcROC{k,1}, lsopcROC{k,2}, lsopcAUC(k), lsopcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPlsPCs, adPlsPCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPlsPCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [plspcROC{k,1}, plspcROC{k,2}, plspcAUC(k), plspcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnWCSs, adWCSs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadWCSsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [wcsROC{k,1}, wcsROC{k,2}, wcsAUC(k), wcsACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnGCs, adGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [gcROC{k,1}, gcROC{k,2}, gcAUC(k), gcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPGCs, adPGCs, 1, 1);        % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPGCsUtP, topNum);                                 % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pgcROC{k,1}, pgcROC{k,2}, pgcAUC(k), pgcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);        % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnTEs, adTEs, 1, 1);        % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadTEsUtP, topNum);                                % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [teROC{k,1}, teROC{k,2}, teAUC(k), teACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLs, adDLs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadDLsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlROC{k,1}, dlROC{k,2}, dlAUC(k), dlACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs, adDLWs, 1, 1);         % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadDLWsUtP, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k), dlwACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCDLs, adPCDLs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPCDLsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pcdlROC{k,1}, pcdlROC{k,2}, pcdlAUC(k), pcdlACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCDLWs, adPCDLWs, 1, 1);         % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPCDLWsUtP, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pcdlwROC{k,1}, pcdlwROC{k,2}, pcdlwAUC(k), pcdlwACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLGs, adDLGs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadDLGsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [dlgROC{k,1}, dlgROC{k,2}, dlgAUC(k), dlgACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCSs, adPCSs, 1, 1);         % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPCSsUtP, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pcsROC{k,1}, pcsROC{k,2}, pcsAUC(k), pcsACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
    %{        
            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnCPCs, adCPCs, 1, 1);         % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadCPCsUtP, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [cpcROC{k,1}, cpcROC{k,2}, cpcAUC(k), cpfcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnFGESs, adFGESs, 1, 1);         % replece cn*s, ad*s
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadFGESsUtP, topNum);                                  % replace cnad*sUtP
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [fgesROC{k,1}, fgesROC{k,2}, fgesAUC(k), fgesACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
    %}
            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnFCAs, adFCAs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadFCAsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [fcaROC{k,1}, fcaROC{k,2}, fcaAUC(k), fcaACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnTSFCs, adTSFCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadTsFCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [tsfcROC{k,1}, tsfcROC{k,2}, tsfcAUC(k), tsfcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnTSFCAs, adTSFCAs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadTsFCAsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [tsfcaROC{k,1}, tsfcaROC{k,2}, tsfcaAUC(k), tsfcaACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMVARDIs, adMVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMvarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mvardiROC{k,1}, mvardiROC{k,2}, mvardiAUC(k), mvardiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPVARDIs, adPVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPvarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pvardiROC{k,1}, pvardiROC{k,2}, pvardiAUC(k), pvardiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMPCVARDIs, adMPCVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMpcvarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mpcvardiROC{k,1}, mpcvardiROC{k,2}, mpcvardiAUC(k), mpcvardiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMPCVARGCs, adMPCVARGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMpcvarGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mpcvargcROC{k,1}, mpcvargcROC{k,2}, mpcvargcAUC(k), mpcvargcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPPCVARDIs, adPPCVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPpcvarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [ppcvardiROC{k,1}, ppcvardiROC{k,2}, ppcvardiAUC(k), ppcvardiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPPCVARGCs, adPPCVARGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPpcvarGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [ppcvargcROC{k,1}, ppcvargcROC{k,2}, ppcvargcAUC(k), ppcvargcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMPLSVARDIs, adMPLSVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMplsvarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mplsdiROC{k,1}, mplsdiROC{k,2}, mplsdiAUC(k), mplsdiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMPLSVARGCs, adMPLSVARGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMplsvarGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mplsgcROC{k,1}, mplsgcROC{k,2}, mplsgcAUC(k), mplsgcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPPLSVARDIs, adPPLSVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPplsvarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pplsdiROC{k,1}, pplsdiROC{k,2}, pplsdiAUC(k), pplsdiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPPLSVARGCs, adPPLSVARGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPplsvarGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pplsgcROC{k,1}, pplsgcROC{k,2}, pplsgcAUC(k), pplsgcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMLSOVARDIs, adMLSOVARDIs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMlsovarDIsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mlsodiROC{k,1}, mlsodiROC{k,2}, mlsodiAUC(k), mlsodiACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnMLSOVARGCs, adMLSOVARGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadMlsovarGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [mlsogcROC{k,1}, mlsogcROC{k,2}, mlsogcAUC(k), mlsogcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);

            i = i + 1;
            [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnPCGCs, adPCGCs, 1, 1);
            [pvList(:,i), I, X] = sortAndPairPValues(control, target, cnadPCGCsUtP, topNum);
            sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
            [pcgcROC{k,1}, pcgcROC{k,2}, pcgcAUC(k), pcgcACC{k}] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);
        end
    end
%    figure; boxplot(X);

    % save result
    fname = [resultsPath '/' resultsPrefix '-cn-ad-roi' num2str(132) '-result.mat'];
    save(fname, 'cosSim', 'fcAUC','pcAUC','pcpcAUC','lsopcAUC','plspcAUC','wcsAUC','gcAUC','pgcAUC','dlAUC','dlwAUC','pcdlAUC','pcdlwAUC','dlgAUC','teAUC','pcsAUC','cpcAUC','fgesAUC','fcaAUC','tsfcAUC','tsfcaAUC','mvardiAUC','mpcvardiAUC','mplsdiAUC','mlsodiAUC','mpcvargcAUC','mplsgcAUC','mlsogcAUC','pcgcAUC', ...
        'fcROC','pcROC','pcpcROC','lsopcROC','plspcROC','wcsROC','gcROC','pgcROC','dlROC','dlwROC','pcdlROC','pcdlwROC','dlgROC','teROC','pcsROC','cpcROC','fgesROC','fcaROC','tsfcROC','tsfcaROC','mvardiROC','mpcvardiROC','mplsdiROC','mlsodiROC','mpcvargcROC','mplsgcROC','mlsogcROC','pcgcROC', ...
        'fcACC','pcACC','pcpcACC','lsopcACC','plspcACC','wcsACC','gcACC','pgcACC','dlACC','dlwACC','pcdlACC','pcdlwACC','dlgACC','teACC','pcsACC','cpcACC','fgesACC','fcaACC','tsfcACC','tsfcaACC','mvardiACC','mpcvardiACC','mplsdiACC','mlsodiACC','mpcvargcACC','mplsgcACC','mlsogcACC','pcgcACC', ...
        'sigCntCN', 'sigCntAD');
    disp('AUCs');
    mean(dlAUC) % show result AUC
    mean(dlwAUC) % show result AUC
    mean(fcAUC) % show result AUC
    mean(pgcAUC) % show result AUC
    mean(pcAUC) % show result AUC
    mean(pcpcAUC) % show result AUC
    mean(wcsAUC) % show result AUC
    mean(dlgAUC) % show result AUC
    mean(pcsAUC) % show result AUC
    mean(fcaAUC) % show result AUC
    mean(tsfcAUC) % show result AUC
    mean(tsfcaAUC) % show result AUC
    mean(mvardiAUC) % show result AUC
    mean(pvardiAUC) % show result AUC
    mean(mpcvardiAUC) % show result AUC
    mean(mpcvargcAUC) % show result AUC
    mean(ppcvardiAUC) % show result AUC
    mean(ppcvargcAUC) % show result AUC
    mean(mlsodiAUC) % show result AUC
    mean(mlsogcAUC) % show result AUC
    mean(pcgcAUC) % show result AUC
    
    % show accuracy
    disp('ACCs');
    fcMaxACC = nan(N,1);
    dlMaxACC = nan(N,1);
    dlwMaxACC = nan(N,1);
    for k=1:N
        fcMaxACC(k) = max(fcACC{k});
        dlMaxACC(k) = max(dlACC{k});
        dlwMaxACC(k) = max(dlwACC{k});
    end
    mean(dlMaxACC) % show result ACC
    mean(dlwMaxACC) % show result ACC
    mean(fcMaxACC) % show result ACC
    
    % plot ROC curve
    figure;
    hold on;
    plotErrorROCcurve(fcROC, N, [0.8,0.2,0.2]);
    plotErrorROCcurve(pcROC, N, [0.8,0.2,0.2]);
%    plotErrorROCcurve(wcsROC, N, [0.9,0.5,0]);
    plotErrorROCcurve(gcROC, N, [0.2,0.8,0.2]);
%    plotErrorROCcurve(pgcROC, N, [0.0,0.5,0.0]);
    plotErrorROCcurve(dlROC, N, [0.2,0.2,0.2]);
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % TODO:
%    plotErrorROCcurve(dcmROC, N, [0.2,0.2,0.8]);
%    plotErrorROCcurve(rnnROC, N, [0.8,0.8,0.2]);
    plotErrorROCcurve(teROC, N, [0.2,0.6,0.8]);
%    plotErrorROCcurve(nnnueROC, N, [0.8,0.2,0.8]);
%    plotErrorROCcurve(dlgROC, N, [0.6,0.6,0.3]);
%    plotErrorROCcurve(pcsROC, N, [0.5,0.5,0.5]);
%    plotErrorROCcurve(cpcROC, N, [0.5,0.5,0.5]);
%    plotErrorROCcurve(fgesROC, N, [0.5,0.5,0.5]);
    plotAverageROCcurve(fcROC, N, '-', [0.8,0.2,0.2],0.5);
    plotAverageROCcurve(pcROC, N, '-', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(pcpcROC, N, '--', [0.5,0.1,0.1],0.5);
    plotAverageROCcurve(lsopcROC, N, '-.', [0.5,0.1,0.1],0.5);
%    plotAverageROCcurve(wcsROC, N, '--', [0.9,0.5,0],0.5);
    plotAverageROCcurve(gcROC, N, '-', [0.1,0.8,0.1],0.5);
%    plotAverageROCcurve(pgcROC, N, '--', [0.0,0.5,0.0],0.5);
%    plotAverageROCcurve(dlROC, N, '-', [0.2,0.2,0.2],1.2);
%    plotAverageROCcurve(dlwROC, N, '--', [0.2,0.2,0.2],0.7); % TODO:
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.2);
    plotAverageROCcurve(dlROC, N, '--', [0.2,0.2,0.2],0.7); % TODO:
%    plotAverageROCcurve(dcmROC, N, '-', [0.2,0.2,0.8],0.5);
%    plotAverageROCcurve(rnnROC, N, '--', [0.7,0.7,0.2],0.5);
    plotAverageROCcurve(teROC, N, '--', [0.2,0.5,0.7],0.5);
%    plotAverageROCcurve(nnnueROC, N, '--', [0.7,0.2,0.7],0.5);
%    plotAverageROCcurve(dlgROC, N, '-.', [0.6,0.6,0.3],0.5);
%    plotAverageROCcurve(pcsROC, N, '-', [0.5,0.5,0.5],0.5);
%    plotAverageROCcurve(cpcROC, N, '--', [0.5,0.5,0.5],0.5);
%    plotAverageROCcurve(fgesROC, N, '-.', [0.5,0.5,0.5],0.5);
%    plotAverageROCcurve(fcaROC, N, '-.', [0.8,0.2,0.2],0.5);
%    plotAverageROCcurve(tsfcROC, N, '-', [0.6,0.2,0.2],1.2);
%    plotAverageROCcurve(tsfcaROC, N, '-.', [0.6,0.2,0.2],1.2);
    plotAverageROCcurve(mvardiROC, N, '-', [0.3,0.3,0.3],0.5);
    plotAverageROCcurve(pvardiROC, N, '--', [0.3,0.3,0.3],0.5);
    plotAverageROCcurve(mpcvardiROC, N, '-', [0.3,0.6,0.6],1.0);
    plotAverageROCcurve(mpcvargcROC, N, '--', [0.3,0.6,0.6],0.8);
    plotAverageROCcurve(ppcvardiROC, N, '-', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(ppcvargcROC, N, '--', [0.3,0.6,0.6],0.5);
    plotAverageROCcurve(mplsdiROC, N, '-', [0.7,0.9,0.9],1.0);
    plotAverageROCcurve(mplsgcROC, N, '--', [0.7,0.9,0.9],0.8);
%    plotAverageROCcurve(pplsdiROC, N, '-', [0.7,0.9,0.9],0.5);
%    plotAverageROCcurve(pplsgcROC, N, '--', [0.7,0.9,0.9],0.5);
    plotAverageROCcurve(mlsodiROC, N, '-.', [0.7,0.9,0.9],1.0);
    plotAverageROCcurve(mlsogcROC, N, ':', [0.7,0.9,0.9],0.8);
    plotAverageROCcurve(pcgcROC, N, '-.', [0.3,0.6,0.6],0.5);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve']);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

function [ampDiffROI,gAmpSigma,pAmpSigma] = checkAmplitudeDiffROI(si1, si2)
    n = size(si1{1},1);
    gAmpSigma = nan(n,length(si1));
    pAmpSigma = nan(n,length(si2));
    for i=1:length(si1)
        gAmpSigma(:,i) = std(si1{i},1,2);
    end
    for i=1:length(si2)
        pAmpSigma(:,i) = std(si2{i},1,2);
    end
    ampDiffROI = nan(n,1);
    for i=1:n
        [ampDiffROI(i), h] = ranksum(gAmpSigma(i,:),pAmpSigma(i,:));
    end
    figure; bar(ampDiffROI); title(['ampDiffROI : si1 vs si2 : P-value']);
%%{
    [B,I]=sort(ampDiffROI);
    for i=1:1
        eg=[0:0.05:2];
        figure; hold on; histogram(gAmpSigma(I(i),:),eg); histogram(pAmpSigma(I(i),:),eg); hold off;
        title(['ampDiff top' num2str(i) ' ROI(' num2str(I(i)) ')']);
    end
%%}
end

function boxplotTwoGroup(X, Y)
    figure; boxplot([X(:), Y(:)]);
end

function [c, r] = index2rowCol(idx,n)
    c=1+floor((idx-1)/n); r=1+mod(idx-1,n);
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

function [x, y, auc, accuracy] = calcAlzROCcurve(control, target, start)
    x = [0]; y = [0]; % start from (0,0)
    accuracy = nan(start+1,1);
    tpmax = length(control);
    fpmax = length(target);
    for i=start:-1:0
        tp = length(find(control>=i));
        fp = length(find(target>=i));
        tn = fpmax - fp;
        x = [x fp/fpmax];
        y = [y tp/tpmax];
        accuracy(i+1) = (tp + tn) / (tpmax + fpmax);
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
%    figure; boxplot(X.');
%    figure; bar(sigCount);
end

function [normalities, normalitiesP] = calculateAlzNormalityTest(weights, roiNames, group, algorithm)
    % constant value
    ROINUM = size(weights,1);

    global resultsPath;
    global resultsPrefix;
    outfName = [resultsPath '/' resultsPrefix '-' algorithm '-' group '-roi' num2str(ROINUM) '-normality.mat'];
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

    global resultsPath;
    global resultsPrefix;
    outfName = [resultsPath '/' resultsPrefix '-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROINUM) '-utest.mat'];
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

    return; % do not show figures
    
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
