% this function works after analyzeAlzheimerDLCM.m

function simulateAlzheimerDLCM2
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');

    cnSbjNum = length(cnSignals);
    adSbjNum = length(adSignals);
    nodeNum = size(cnSignals{1},1);
    
    % get DLCM-GC from original CN and AD (normal and recovery training)
    [cnDLs, meanCnDL, stdCnDL] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcm', 1);
    [adDLs, meanAdDL, stdAdDL] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcm', 1);
    sigCnDLs = (cnDLs - nanmean(cnDLs(:))) / nanstd(cnDLs(:),1);
    sigAdDLs = (adDLs - nanmean(adDLs(:))) / nanstd(adDLs(:),1);

    [rccnDLs, meanRcCnDL, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'dlcmrc', 1); % do recovery training
    [rcadDLs, meanRcAdDL, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'dlcmrc', 1); % do recovery training

    % simulate CN & AD signals from first frame
    [cnDLWs, smcnSignals, cnSubDLWs] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlw', 'cn');
    [adDLWs, smadSignals, adSubDLWs] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlw', 'ad');
    sigCnDLWs = (cnDLWs - nanmean(cnDLWs(:))) / nanstd(cnDLWs(:),1);
    sigAdDLWs = (adDLWs - nanmean(adDLWs(:))) / nanstd(adDLWs(:),1);
    meanCnDLW = nanmean(cnDLWs,3);
    meanAdDLW = nanmean(adDLWs,3);

    % simulate recovery trained CN & AD signals from first frame
    [rccnDLWs, smrccnSignals] = simulateNodeSignals(cnSignals, roiNames, 'cn', 'dlwrc', 'cn');
    [rcadDLWs, smrcadSignals] = simulateNodeSignals(adSignals, roiNames, 'ad', 'dlwrc', 'ad');

    % simulate CN & AD signals from last-1 frame
    cn2Signals = cnSignals;
    ad2Signals = adSignals;
    for i=1:cnSbjNum
        cn2Signals{i}(:,1) = cn2Signals{i}(:,end-1);
        cn2Signals{i}(:,2) = cn2Signals{i}(:,end);
    end
    for i=1:adSbjNum
        ad2Signals{i}(:,1) = ad2Signals{i}(:,end-1);
        ad2Signals{i}(:,2) = ad2Signals{i}(:,end);
    end
    [~, smcn2Signals] = simulateNodeSignals(cn2Signals, roiNames, 'cn2', 'dlw', 'cn');
    [~, smad2Signals] = simulateNodeSignals(ad2Signals, roiNames, 'ad2', 'dlw', 'ad');

    % expanding amplitude of simulated CN & AD signals from last-1 frame (type1)
    % -- no effect
%{
    smcn3Signals = smcn2Signals;
    smad3Signals = smad2Signals;
    for i=1:cnSbjNum
        smcn3Signals{i} = expandAmplitude(smcn2Signals{i}, 1.3);
    end
    for i=1:adSbjNum
        smad3Signals{i} = expandAmplitude(smad2Signals{i}, 1.3);
    end
%}
    % expanding amplitude of simulated CN & AD signals from last-1 frame (type2)
    % finding best rate -- result showed bigger k is better (k=10 is best)
%{
    smcn4Signals = smcn2Signals;
    smad4Signals = smad2Signals;
    cnsmcn4DLWrs = zeros(1,10);
    for k=10:10
        for i=1:cnSbjNum
            smcn4Signals{i} = expandAmplitude2(smcn2Signals{i}, 1+0.5*k);
        end
        for i=1:adSbjNum
            smad4Signals{i} = expandAmplitude2(smad2Signals{i}, 1+0.5*k);
        end
        name = ['smcn4_' num2str(k)];
        [smcn4DLs, meanSmcn4DL, ~] = calculateConnectivity(smcn4Signals, roiNames, name, 'dlcm', 1);
        [smcn4DLWs, meanSmcn4DLW, ~, smcn4SubDLWs] = calculateConnectivity(smcn4Signals, roiNames, name, 'dlw', 1);
        figure; cnsmcn4DLWrs(k) = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn4DLW);
    end
%}
    % re-train CN signals with expanded EC amplitude (type3)
    % -- no effect. trained DLCM network showed bad simulation signals
    S2 = ones(nodeNum, nodeNum+1);
    S2(:,2:end) = S2(:, 2:end) - eye(nodeNum);
    IS2 = ones(nodeNum, nodeNum+1);
    %[smcn5DLWs, smcn5Signals] = retrainDLCMAndECmultiPattern(cnSignals, cnDLWs, cnSubDLWs, S2, IS2, roiNames, 'cn5');

    % --------------------------------------------------------------------------------------------------------------
    % check DLCM-EC and DLCM-GC of simulated CN and AD
    [smcnDLs, meanSmcnDL, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlcm', 1);
    [smcnDLWs, meanSmcnDLW, ~, smcnSubDLWs] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'dlw', 1);
    [smadDLs, meanSmadDL, ~] = calculateConnectivity(smadSignals, roiNames, 'smad', 'dlcm', 1);
    [smadDLWs, meanSmadDLW, ~, smadSubDLWs] = calculateConnectivity(smadSignals, roiNames, 'smad', 'dlw', 1);
    sigSmcnDLWs = (smcnDLWs - nanmean(smcnDLWs(:))) / nanstd(smcnDLWs(:),1);
    sigSmadDLWs = (smadDLWs - nanmean(smadDLWs(:))) / nanstd(smadDLWs(:),1);

    [smrccnDLs, meanSmrccnDL, ~] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'dlcm', 1);
    [smrccnDLWs, meanSmrccnDLW, ~, smrccnSubDLWs] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'dlw', 1);
    [smrcadDLs, meanSmrcadDL, ~] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'dlcm', 1);
    [smrcadDLWs, meanSmrcadDLW, ~, smrcadSubDLWs] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'dlw', 1);
%{
    for i=1:cnSbjNum
        plotCorrelationZiZij([], cnSubDLWs(:,:,i), [], smcnSubDLWs(:,:,i), nodeNum, ['sbj' num2str(i)], 'original', 'smcn');
        plotCorrelationZiZij([], cnSubDLWs(:,:,i), [], smrccnSubDLWs(:,:,i), nodeNum, ['sbj' num2str(i)], 'original', 'smrccn');
    end
%}
    [smcn2DLs, meanSmcn2DL, ~] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'dlcm', 1);
    [smcn2DLWs, meanSmcn2DLW, ~, smcn2SubDLWs] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'dlw', 1);
    [smad2DLs, meanSmad2DL, ~] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'dlcm', 1);
    [smad2DLWs, meanSmad2DLW, ~, smad2SubDLWs] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'dlw', 1);
    sigSmcn2DLWs = (smcn2DLWs - nanmean(smcn2DLWs(:))) / nanstd(smcn2DLWs(:),1);
    sigSmad2DLWs = (smad2DLWs - nanmean(smad2DLWs(:))) / nanstd(smad2DLWs(:),1);

%    [smcn3DLs, meanSmcn3DL, ~] = calculateConnectivity(smcn3Signals, roiNames, 'smcn3', 'dlcm', 1);
%    [smcn3DLWs, meanSmcn3DLW, ~] = calculateConnectivity(smcn3Signals, roiNames, 'smcn3', 'dlw', 1);
%    [smad3DLs, meanSmad3DL, ~] = calculateConnectivity(smad3Signals, roiNames, 'smad3', 'dlcm', 1);
%    [smad3DLWs, meanSmad3DLW, ~] = calculateConnectivity(smad3Signals, roiNames, 'smad3', 'dlw', 1);

%    [smcn4DLs, meanSmcn4DL, ~] = calculateConnectivity(smcn4Signals, roiNames, 'smcn4', 'dlcm', 1);
%    [smcn4DLWs, meanSmcn4DLW, ~] = calculateConnectivity(smcn4Signals, roiNames, 'smcn4', 'dlw', 1);
%    [smad4DLs, meanSmad4DL, ~] = calculateConnectivity(smad4Signals, roiNames, 'smad4', 'dlcm', 1);
%    [smad4DLWs, meanSmad4DLW, ~] = calculateConnectivity(smad4Signals, roiNames, 'smad4', 'dlw', 1);

    % check relation between Zi vs signal mean diff, and Zij vs signal amplitude (change teaching signal)
%    checkRelationSubDLWandSignals(cnSignals, cnDLWs, cnSubDLWs, 'cn', 0);
%    checkRelationSubDLWandSignals(smcnSignals, smcnDLWs, smcnSubDLWs, 'smcn', 1);
%    checkRelationSubDLWandSignals2(cnSignals, cnDLWs, cnSubDLWs, smcnDLWs, smcnSubDLWs, 'cn');

    % -- move mean (range=65) based amplitude expansion of whole signal
    % -- change original cnSignals -> simulate -> calc EC
    % -- somehow working well
    [smcn7DLWs, smcn7SubDLWs, smcn7Signals, smcn7DLs] = checkRelationSubDLWandSignals3(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 1);
    [smcn8DLWs, smcn8SubDLWs, smcn8Signals, smcn8DLs] = checkRelationSubDLWandSignals3b(cn2Signals, cnDLWs, cnSubDLWs, smcn2Signals, smcn2DLWs, smcn2SubDLWs, smcn2DLs, 'cn2', 'cn');
    meanSmcn7DLW = nanmean(smcn7DLWs,3);
    meanSmcn8DLW = nanmean(smcn8DLWs,3);
    meanSmcn7DL = nanmean(smcn7DLs,3);
    meanSmcn8DL = nanmean(smcn8DLs,3);
    sigSmcn7DLWs = (smcn7DLWs - nanmean(smcn7DLWs(:))) / nanstd(smcn7DLWs(:),1);
    sigSmcn8DLWs = (smcn8DLWs - nanmean(smcn8DLWs(:))) / nanstd(smcn8DLWs(:),1);
    [smad7DLWs, smad7SubDLWs, smad7Signals, smad7DLs] = checkRelationSubDLWandSignals3(adSignals, adDLWs, adSubDLWs, smadSignals, smadDLWs, smadSubDLWs, smadDLs, 'ad', 1);
    [smad8DLWs, smad8SubDLWs, smad8Signals, smad8DLs] = checkRelationSubDLWandSignals3b(ad2Signals, adDLWs, adSubDLWs, smad2Signals, smad2DLWs, smad2SubDLWs, smad2DLs, 'ad2', 'ad');
    meanSmad7DLW = nanmean(smad7DLWs,3);
    meanSmad8DLW = nanmean(smad8DLWs,3);
    meanSmad7DL = nanmean(smad7DLs,3);
    meanSmad8DL = nanmean(smad8DLs,3);
    sigSmad7DLWs = (smad7DLWs - nanmean(smad7DLWs(:))) / nanstd(smad7DLWs(:),1);
    sigSmad8DLWs = (smad8DLWs - nanmean(smad8DLWs(:))) / nanstd(smad8DLWs(:),1);
    % -- type=2 does not perform better than type=1
%    [smcn11DLWs, smcn11SubDLWs, smcn11Signals, smcn11DLs] = checkRelationSubDLWandSignals3(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 2);
    % -- type=3 does not perform better than type=1
%    [smcn12DLWs, smcn12SubDLWs, smcn12Signals, smcn12DLs] = checkRelationSubDLWandSignals3(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 3);
    % -- type=4 extend amplitudes of high frequency (1:15) via wavelet transfrom & invert wavelet transfrom
    % -- this does not perform better than type=1
%    [smcn13DLWs, smcn13SubDLWs, smcn13Signals, smcn13DLs] = checkRelationSubDLWandSignals3(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 4);
    % -- type=5 extend amplitudes of low frequency via wavelet transfrom & invert wavelet transfrom
    % -- this does not work 
%    [smcn14DLWs, smcn14SubDLWs, smcn14Signals, smcn14DLs] = checkRelationSubDLWandSignals3(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 5);
    % -- type=6 extend amplitudes of high frequency (1:3) via wavelet transfrom & invert wavelet transfrom
    % -- this does not perform better than type=1 and type=4
%    [smcn15DLWs, smcn15SubDLWs, smcn15Signals, smcn15DLs] = checkRelationSubDLWandSignals3(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 6);

    % -- check wavelet transform effect
    % -- change original cnSignals -(cwt)-> frequency -(icwt+movemean)-> Signals -> calc EC
    [wtcnDLWs, wtcnSubDLWs, wtcnSignals, wtcnDLs] = checkWaveletTransformEffect(cnSignals, cnDLWs, cnSubDLWs, 'cn', 1);
    meanWtcnDLW = nanmean(wtcnDLWs,3);
    meanWtcnDL = nanmean(wtcnDLs,3);
    % -- change original cnSignals -(stft)-> frequency -(istft)-> Signals -> calc EC
%    [wtcn2DLWs, wtcn2SubDLWs, wtcn2Signals, wtcn2DLs] = checkWaveletTransformEffect(cnSignals, cnDLWs, cnSubDLWs, 'cn', 2);
%    meanWtcn2DLW = nanmean(wtcn2DLWs,3);
%    meanWtcn2DL = nanmean(wtcn2DLs,3);
    % -- change original cnSignals -(cwt)-> frequency(remove high) -(icwt+movemean)-> Signals -> calc EC
    [wtcn3DLWs, wtcn3SubDLWs, wtcn3Signals, wtcn3DLs] = checkWaveletTransformEffect(cnSignals, cnDLWs, cnSubDLWs, 'cn', 3);
    meanWtcn3DLW = nanmean(wtcn3DLWs,3);
    meanWtcn3DL = nanmean(wtcn3DLs,3);
    % -- change original cnSignals -(cwt)-> frequency -(icwt)-> Signals -> calc EC
    [wtcn4DLWs, wtcn4SubDLWs, wtcn4Signals, wtcn4DLs] = checkWaveletTransformEffect(cnSignals, cnDLWs, cnSubDLWs, 'cn', 4);
    meanWtcn4DLW = nanmean(wtcn4DLWs,3);
    meanWtcn4DL = nanmean(wtcn4DLs,3);

    % check relation between Zij vs signal amplitude (change other input signals)
    % -- not working well
%    checkRelationSubDLWandSignals4(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, 'smcn');

    % check relation between Zij vs node weights (change other input signals)
    % -- very effective, but only for network weight
%    checkRelationSubDLWandWeights(cnSignals, cnSubDLWs, smcnSubDLWs, 'smcn', 1);
%    checkRelationSubDLWandWeights(cnSignals, cnSubDLWs, smcnSubDLWs, 'smcn', 2);

    % check relation between Zij vs node weights (change other input signals)
    % -- change original DLCM weight -> simulate -> calc EC
    % -- somehow working well, but not better than smcn7DLWs
    [smcn9DLWs, smcn9SubDLWs, smcn9Signals, smcn9DLs] = checkRelationSubDLWandWeights2(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 1);
    meanSmcn9DLW = nanmean(smcn9DLWs,3);
    meanSmcn9DL = nanmean(smcn9DLs,3);
    % -- type=2 does not work at all
%    [smcn10DLWs, smcn10SubDLWs, smcn10Signals, smcn10DLs] = checkRelationSubDLWandWeights2(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, smcnDLs, 'cn', 2);
%    meanSmcn10DLW = nanmean(smcn10DLWs,3);
%    meanSmcn10DL = nanmean(smcn10DLs,3);

    % re-train CN signals with shifting signals and expanding EC amplitude (type4)
    % -- parallel shift does not affect EC (both Zi and Zij shifted)
    % -- move mean (range=5) based amplitude expansion does not work well
    % -- change simulated signals -> calc EC
    [smcn6DLWs, smcn6SubDLWs, smcn6Signals] = shiftAndExpandAmplitude(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, 'smcn', 1);
    meanSmcn6DLW = nanmean(smcn6DLWs,3);
    % -- change simulated signals of high frequency (1:3) (wavelet) -> calc EC
    % -- this works a bit.
    [smcn6bDLWs, smcn6bSubDLWs, smcn6bSignals] = shiftAndExpandAmplitude(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, 'smcn', 3);
    meanSmcn6bDLW = nanmean(smcn6bDLWs,3);
    % -- change simulated signals of high frequency (1:10) (wavelet) -> calc EC
    % -- 
%    [smcn6cDLWs, smcn6cSubDLWs, smcn6cSignals] = shiftAndExpandAmplitude(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, 'smcn', 4);
%    meanSmcn6cDLW = nanmean(smcn6cDLWs,3);
    % -- 
    [smcn6dDLWs, smcn6dSubDLWs, smcn6dSignals] = shiftAndExpandAmplitude(cnSignals, cnDLWs, cnSubDLWs, smcnSignals, smcnDLWs, smcnSubDLWs, 'smcn', 5);

%{
    figure; cnsmcnDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcnDLW);
    figure; adsmadDLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmadDLW);
    figure; cnsmrccnDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmrccnDLW);
    figure; adsmrcadDLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmrcadDLW);
    figure; cnsmcn2DLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn2DLW);
    figure; adsmad2DLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmad2DLW);
    
    m1 = nanmean(meanCnDLW(:));
    m2 = nanmean(meanSmcnDLW(:));
    s1 = nanstd(meanCnDLW(:));
    s2 = nanstd(meanSmcnDLW(:));
    sigCnDLW = (meanCnDLW - m1) ./ s1;
    sigSmcnDLW = (meanSmcnDLW - m2) ./ s2;
    figure; cnsmcnDLWr = plotTwoSignalsCorrelation(sigCnDLW, sigSmcnDLW);
%}    
%    figure; cnsmcn3DLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn3DLW);
%    figure; adsmad3DLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmad3DLW);
%    figure; cnsmcn4DLWr = plotTwoSignalsCorrelation(meanCnDLW, meanSmcn4DLW);
%    figure; adsmad4DLWr = plotTwoSignalsCorrelation(meanAdDLW, meanSmad4DLW);

    % check FC of simulated CN and AD
    [cnFCs, meanCnFC, ~] = calculateConnectivity(cnSignals, roiNames, 'cn', 'fc', 1);
    [adFCs, meanAdFC, ~] = calculateConnectivity(adSignals, roiNames, 'ad', 'fc', 1);
    [smcnFCs, meanSmcnFC, ~] = calculateConnectivity(smcnSignals, roiNames, 'smcn', 'fc', 1);
    [smadFCs, meanSmadFC, ~] = calculateConnectivity(smadSignals, roiNames, 'smad', 'fc', 1);
%    [smrccnFCs, meanSmrccnFC, ~] = calculateConnectivity(smrccnSignals, roiNames, 'smrccn', 'fc', 1);
%    [smrcadFCs, meanSmrcadFC, ~] = calculateConnectivity(smrcadSignals, roiNames, 'smrcad', 'fc', 1);
    [smcn2FCs, meanSmcn2FC, ~] = calculateConnectivity(smcn2Signals, roiNames, 'smcn2', 'fc', 1);
    [smad2FCs, meanSmad2FC, ~] = calculateConnectivity(smad2Signals, roiNames, 'smad2', 'fc', 1);
    [smcn6FCs, meanSmcn6FC, ~] = calculateConnectivity(smcn6Signals, roiNames, 'smcn6', 'fc', 1);
    [smcn7FCs, meanSmcn7FC, ~] = calculateConnectivity(smcn7Signals, roiNames, 'smcn7', 'fc', 1);
    [smcn8FCs, meanSmcn8FC, ~] = calculateConnectivity(smcn8Signals, roiNames, 'smcn8', 'fc', 1);
    [smcn9FCs, meanSmcn9FC, ~] = calculateConnectivity(smcn9Signals, roiNames, 'smcn9', 'fc', 1);
    [smad7FCs, meanSmad7FC, ~] = calculateConnectivity(smad7Signals, roiNames, 'smad7', 'fc', 1);
    [smad8FCs, meanSmad8FC, ~] = calculateConnectivity(smad8Signals, roiNames, 'smad8', 'fc', 1);
    [wtcnFCs, meanWtcnFC, ~] = calculateConnectivity(wtcnSignals, roiNames, 'wtcn', 'fc', 1);
%    [wtcn2FCs, meanWtcn2FC, ~] = calculateConnectivity(wtcn2Signals, roiNames, 'wtcn2', 'fc', 1);
    [wtcn3FCs, meanWtcn3FC, ~] = calculateConnectivity(wtcn3Signals, roiNames, 'wtcn3', 'fc', 1);
    [wtcn4FCs, meanWtcn4FC, ~] = calculateConnectivity(wtcn4Signals, roiNames, 'wtcn4', 'fc', 1);
    
%{
    figure; cnsmcnFCr = plotTwoSignalsCorrelation(meanCnFC, meanSmcnFC);
    figure; adsmadFCr = plotTwoSignalsCorrelation(meanAdFC, meanSmadFC);
    figure; cnsmrccnFCr = plotTwoSignalsCorrelation(meanCnFC, meanSmrccnFC);
    figure; adsmrcadFCr = plotTwoSignalsCorrelation(meanAdFC, meanSmrcadFC);
    figure; cnsmcn2FCr = plotTwoSignalsCorrelation(meanCnFC, meanSmcn2FC);
    figure; adsmad2FCr = plotTwoSignalsCorrelation(meanAdFC, meanSmad2FC);
%}

    % --------------------------------------------------------------------------------------------------------------
    % PCA test -- not work well. first 3 components describe only 15% of whole variable.
    % -- and first 3 most significantly different HC vs AD components are still not separated well.
    [score, explained, pvals] = checkPCAHCvsAD(cnDLWs, adDLWs, nodeNum, 'DLCM-EC');
    [score, explained, pvals] = checkPCAHCvsAD(cnDLs, adDLs, nodeNum, 'DLCM-GC');
    [score, explained, pvals] = checkPCAHCvsAD(cnFCs, adFCs, nodeNum, 'FC');
    [score, explained, pvals] = checkPCAHCvsAD(smcnDLWs, smadDLWs, nodeNum, 'DLCM-EC');
    
    % --------------------------------------------------------------------------------------------------------------
    % plot correlation and cos similarity
    algNum = 30;
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    cosSim(2) = getCosSimilarity(meanCnDLW, meanSmcn2DLW);
    cosSim(3) = getCosSimilarity(meanCnDLW, meanSmcn6DLW);
    cosSim(4) = getCosSimilarity(meanCnDLW, meanSmcn7DLW);
    cosSim(5) = getCosSimilarity(meanCnDLW, meanSmcn8DLW);
    cosSim(6) = getCosSimilarity(meanCnDLW, meanSmcn9DLW);
    cosSim(7) = getCosSimilarity(meanCnDLW, meanWtcnDLW);
%    cosSim(8) = getCosSimilarity(meanCnDLW, meanWtcn2DLW);
    cosSim(8) = getCosSimilarity(meanCnDLW, meanWtcn3DLW);
    cosSim(9) = getCosSimilarity(meanCnDLW, meanWtcn4DLW);
    cosSim(11) = getCosSimilarity(meanCnDL, meanSmcnDL);
    cosSim(12) = getCosSimilarity(meanCnDL, meanSmcn2DL);
    cosSim(13) = getCosSimilarity(meanCnDL, meanSmcn7DL);
    cosSim(14) = getCosSimilarity(meanCnDL, meanSmcn8DL);
    cosSim(15) = getCosSimilarity(meanCnDL, meanSmcn9DL);
    cosSim(16) = getCosSimilarity(meanCnDL, meanWtcnDL);
%    cosSim(17) = getCosSimilarity(meanCnDL, meanWtcn2DL);
    cosSim(17) = getCosSimilarity(meanCnDL, meanWtcn3DL);
    cosSim(18) = getCosSimilarity(meanCnDL, meanWtcn4DL);
    cosSim(21) = getCosSimilarity(meanCnFC, meanSmcnFC);
    cosSim(22) = getCosSimilarity(meanCnFC, meanSmcn2FC);
    cosSim(23) = getCosSimilarity(meanCnFC, meanSmcn6FC);
    cosSim(24) = getCosSimilarity(meanCnFC, meanSmcn7FC);
    cosSim(25) = getCosSimilarity(meanCnFC, meanSmcn8FC);
    cosSim(26) = getCosSimilarity(meanCnFC, meanSmcn9FC);
    cosSim(27) = getCosSimilarity(meanCnFC, meanWtcnFC);
%    cosSim(28) = getCosSimilarity(meanCnFC, meanWtcn2FC);
    cosSim(28) = getCosSimilarity(meanCnFC, meanWtcn3FC);
    cosSim(29) = getCosSimilarity(meanCnFC, meanWtcn4FC);
%    X = categorical({'dlec-cn-smcn','dlec-cn-smcn2','dlec-cn-smcn6','dlec-cn-smcn7',...
%        'dlgc-cn-smcn','dlgc-cn-smcn2','dlgc-cn-smcn7',...
%        'fc-cn-smcn','fc-cn-smcn2','fc-cn-smcn6','fc-cn-smcn7'});
    figure; bar(cosSim); title('cos similarity between mean CN matrix and SimCN by each algorithm');

    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanAdDLW, meanSmadDLW);
    cosSim(2) = getCosSimilarity(meanAdDLW, meanSmad2DLW);
    cosSim(3) = getCosSimilarity(meanAdDLW, meanSmad7DLW);
    cosSim(4) = getCosSimilarity(meanAdDLW, meanSmad8DLW);
    cosSim(11) = getCosSimilarity(meanAdDL, meanSmadDL);
    cosSim(12) = getCosSimilarity(meanAdDL, meanSmad2DL);
    cosSim(13) = getCosSimilarity(meanAdDL, meanSmad7DL);
    cosSim(14) = getCosSimilarity(meanAdDL, meanSmad8DL);
    cosSim(21) = getCosSimilarity(meanAdFC, meanSmadFC);
    cosSim(22) = getCosSimilarity(meanAdFC, meanSmad2FC);
    cosSim(23) = getCosSimilarity(meanAdFC, meanSmad7FC);
    cosSim(24) = getCosSimilarity(meanAdFC, meanSmad8FC);
%    X = categorical({'dlec-ad-smad','dlec-ad-smad2','dlgc-ad-smad','dlgc-ad-smad2','fc-ad-smad','fc-ad-smad2','',''});
    figure; bar(cosSim); title('cos similarity between mean AD matrix and SimAD by each algorithm');
%{
    % plot correlation and cos similarity (all)
    algNum = 24;
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanSmcnDLW);
    cosSim(2) = getCosSimilarity(meanAdDLW, meanSmadDLW);
    cosSim(3) = getCosSimilarity(meanAdDLW, meanSmcnDLW);
    cosSim(4) = getCosSimilarity(meanCnDLW, meanSmadDLW);
    cosSim(5) = getCosSimilarity(meanCnDL, meanSmcnDL);
    cosSim(6) = getCosSimilarity(meanAdDL, meanSmadDL);
    cosSim(7) = getCosSimilarity(meanAdDL, meanSmcnDL);
    cosSim(8) = getCosSimilarity(meanCnDL, meanSmadDL);
    cosSim(9) = getCosSimilarity(meanCnFC, meanSmcnFC);
    cosSim(10) = getCosSimilarity(meanAdFC, meanSmadFC);
    cosSim(11) = getCosSimilarity(meanAdFC, meanSmcnFC);
    cosSim(12) = getCosSimilarity(meanCnFC, meanSmadFC);
%    cosSim(13) = getCosSimilarity(meanCnDLW, meanSmrccnDLW);
%    cosSim(14) = getCosSimilarity(meanAdDLW, meanSmrcadDLW);
%    cosSim(15) = getCosSimilarity(meanAdDLW, meanSmrccnDLW);
%    cosSim(16) = getCosSimilarity(meanCnDLW, meanSmrcadDLW);
%    cosSim(17) = getCosSimilarity(meanCnFC, meanSmrccnFC);
%    cosSim(18) = getCosSimilarity(meanAdFC, meanSmrcadFC);
%    cosSim(19) = getCosSimilarity(meanAdFC, meanSmrccnFC);
%    cosSim(20) = getCosSimilarity(meanCnFC, meanSmrcadFC);
    cosSim(21) = getCosSimilarity(meanCnDLW, meanSmcn2DLW);
    cosSim(22) = getCosSimilarity(meanAdDLW, meanSmad2DLW);
    cosSim(23) = getCosSimilarity(meanAdDLW, meanSmcn2DLW);
    cosSim(24) = getCosSimilarity(meanCnDLW, meanSmad2DLW);
    X = categorical({'cn-smcn','ad-smad','ad-smcn','cn-smad','cn-smcn-dlgc','ad-smad-dlgc','ad-smcn-dlgc','cn-smad-dlgc',...
        'cn-smcn-fc','ad-smad-fc','ad-smcn-fc','cn-smad-fc','cn-smrccn','ad-smrcad','ad-smrccn','cn-smrcad',...
        'cn-smrccn-fc','ad-smrcad-fc','ad-smrccn-fc','cn-smrcad-fc','cn-smcn2','ad-smad2','ad-smcn2','cn-smad2',});
    figure; bar(X, cosSim);
    title('cos similarity between CN and SimCN by each algorithm');
%}
    % plot box-and-whisker plot of cos similarity between mean ec matrix and each subject ec
    cosSims = nan(cnSbjNum,30);
    for i=1:cnSbjNum
        cosSims(i,1) = getCosSimilarity(meanCnDLW, cnDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(meanCnDLW, smcnDLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(meanCnDLW, smcn2DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(meanCnDLW, smcn6DLWs(:,:,i));
        cosSims(i,5) = getCosSimilarity(meanCnDLW, smcn7DLWs(:,:,i));
        cosSims(i,6) = getCosSimilarity(meanCnDLW, smcn8DLWs(:,:,i));
        cosSims(i,7) = getCosSimilarity(meanCnDLW, smcn9DLWs(:,:,i));
        cosSims(i,8) = getCosSimilarity(meanCnDLW, wtcnDLWs(:,:,i));
%        cosSims(i,9) = getCosSimilarity(meanCnDLW, wtcn2DLWs(:,:,i));
        cosSims(i,9) = getCosSimilarity(meanCnDLW, wtcn3DLWs(:,:,i));
        cosSims(i,10) = getCosSimilarity(meanCnDLW, wtcn4DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(meanCnDL, cnDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(meanCnDL, smcnDLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(meanCnDL, smcn2DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(meanCnDL, smcn7DLs(:,:,i));
        cosSims(i,15) = getCosSimilarity(meanCnDL, smcn8DLs(:,:,i));
        cosSims(i,16) = getCosSimilarity(meanCnDL, smcn9DLs(:,:,i));
        cosSims(i,17) = getCosSimilarity(meanCnDL, wtcnDLs(:,:,i));
%        cosSims(i,18) = getCosSimilarity(meanCnDL, wtcn2DLs(:,:,i));
        cosSims(i,18) = getCosSimilarity(meanCnDL, wtcn3DLs(:,:,i));
        cosSims(i,19) = getCosSimilarity(meanCnDL, wtcn4DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(meanCnFC, cnFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(meanCnFC, smcnFCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(meanCnFC, smcn2FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(meanCnFC, smcn6FCs(:,:,i));
        cosSims(i,25) = getCosSimilarity(meanCnFC, smcn7FCs(:,:,i));
        cosSims(i,26) = getCosSimilarity(meanCnFC, smcn8FCs(:,:,i));
        cosSims(i,27) = getCosSimilarity(meanCnFC, smcn9FCs(:,:,i));
        cosSims(i,28) = getCosSimilarity(meanCnFC, wtcnFCs(:,:,i));
%        cosSims(i,29) = getCosSimilarity(meanCnFC, wtcn2FCs(:,:,i));
        cosSims(i,29) = getCosSimilarity(meanCnFC, wtcn3FCs(:,:,i));
        cosSims(i,30) = getCosSimilarity(meanCnFC, wtcn4FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:30));
    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
    [sortCosCnDLW,idxCosCnDLW] = sort(cosSims(:,1),'descend');
    [sortCosCnFC,idxCosCnFC] = sort(cosSims(:,21),'descend');

    cosSims = nan(adSbjNum,30);
    for i=1:adSbjNum
        cosSims(i,1) = getCosSimilarity(meanAdDLW, adDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(meanAdDLW, smadDLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(meanAdDLW, smad2DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(meanAdDLW, smad7DLWs(:,:,i));
        cosSims(i,5) = getCosSimilarity(meanAdDLW, smad8DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(meanAdDL, adDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(meanAdDL, smadDLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(meanAdDL, smad2DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(meanAdDL, smad7DLs(:,:,i));
        cosSims(i,15) = getCosSimilarity(meanAdDL, smad8DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(meanAdFC, adFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(meanAdFC, smadFCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(meanAdFC, smad2FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(meanAdFC, smad7FCs(:,:,i));
        cosSims(i,25) = getCosSimilarity(meanAdFC, smad8FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    [anovaP,tbl,stats] = kruskalwallis(cosSims(:,1:30));
    c = multcompare(stats);
%    [p1,h1] = ranksum(cosSims(:,1),cosSims(:,2));
%    [p2,h2] = ranksum(cosSims(:,1),cosSims(:,3));
    [sortCosAdDLW,idxCosAdDLW] = sort(cosSims(:,1),'descend');
    [sortCosAdFC,idxCosAdFC] = sort(cosSims(:,21),'descend');

    % plot correlation between original EC matrix and simulated signals EC matrix
    for i=1:1 %cnSbjNum
        % calc & plot correlation of original Zi vs simulating Zi
        [corrZi, corrZij] = calcCorrelationZiZij(cnSubDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum);
        plotCorrelationZiZij(cnDLWs(:,:,i), cnSubDLWs(:,:,i), smcnDLWs(:,:,i), smcnSubDLWs(:,:,i), nodeNum, ['sbj' num2str(i)], 'original', 'simulating');
    end

    % plot box-and-whisker plot of cos similarity between original ec matrix and simulated signals ec matrix
    cosSims = nan(cnSbjNum,30);
    for i=1:cnSbjNum
        cosSims(i,1) = getCosSimilarity(cnDLWs(:,:,i), smcnDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(cnDLWs(:,:,i), smcn2DLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(cnDLWs(:,:,i), smcn6DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(cnDLWs(:,:,i), smcn7DLWs(:,:,i));
        cosSims(i,5) = getCosSimilarity(cnDLWs(:,:,i), smcn8DLWs(:,:,i));
        cosSims(i,6) = getCosSimilarity(cnDLWs(:,:,i), smcn9DLWs(:,:,i));
        cosSims(i,7) = getCosSimilarity(cnDLWs(:,:,i), wtcnDLWs(:,:,i));
%        cosSims(i,8) = getCosSimilarity(cnDLWs(:,:,i), wtcn2DLWs(:,:,i));
        cosSims(i,8) = getCosSimilarity(cnDLWs(:,:,i), wtcn3DLWs(:,:,i));
        cosSims(i,9) = getCosSimilarity(cnDLWs(:,:,i), wtcn4DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(cnDLs(:,:,i), smcnDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(cnDLs(:,:,i), smcn2DLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(cnDLs(:,:,i), smcn7DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(cnDLs(:,:,i), smcn8DLs(:,:,i));
        cosSims(i,15) = getCosSimilarity(cnDLs(:,:,i), smcn9DLs(:,:,i));
        cosSims(i,16) = getCosSimilarity(cnDLs(:,:,i), wtcnDLs(:,:,i));
%        cosSims(i,17) = getCosSimilarity(cnDLs(:,:,i), wtcn2DLs(:,:,i));
        cosSims(i,17) = getCosSimilarity(cnDLs(:,:,i), wtcn3DLs(:,:,i));
        cosSims(i,18) = getCosSimilarity(cnDLs(:,:,i), wtcn4DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(cnFCs(:,:,i), smcnFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(cnFCs(:,:,i), smcn2FCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(cnFCs(:,:,i), smcn6FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(cnFCs(:,:,i), smcn7FCs(:,:,i));
        cosSims(i,25) = getCosSimilarity(cnFCs(:,:,i), smcn8FCs(:,:,i));
        cosSims(i,26) = getCosSimilarity(cnFCs(:,:,i), smcn9FCs(:,:,i));
        cosSims(i,27) = getCosSimilarity(cnFCs(:,:,i), wtcnFCs(:,:,i));
%        cosSims(i,28) = getCosSimilarity(cnFCs(:,:,i), wtcn2FCs(:,:,i));
        cosSims(i,28) = getCosSimilarity(cnFCs(:,:,i), wtcn3FCs(:,:,i));
        cosSims(i,29) = getCosSimilarity(cnFCs(:,:,i), wtcn4FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    cosSims = nan(cnSbjNum,30);
    for i=1:adSbjNum
        cosSims(i,1) = getCosSimilarity(adDLWs(:,:,i), smadDLWs(:,:,i));
        cosSims(i,2) = getCosSimilarity(adDLWs(:,:,i), smad2DLWs(:,:,i));
        cosSims(i,3) = getCosSimilarity(adDLWs(:,:,i), smad7DLWs(:,:,i));
        cosSims(i,4) = getCosSimilarity(adDLWs(:,:,i), smad8DLWs(:,:,i));
        cosSims(i,11) = getCosSimilarity(adDLs(:,:,i), smadDLs(:,:,i));
        cosSims(i,12) = getCosSimilarity(adDLs(:,:,i), smad2DLs(:,:,i));
        cosSims(i,13) = getCosSimilarity(adDLs(:,:,i), smad7DLs(:,:,i));
        cosSims(i,14) = getCosSimilarity(adDLs(:,:,i), smad8DLs(:,:,i));
        cosSims(i,21) = getCosSimilarity(adFCs(:,:,i), smadFCs(:,:,i));
        cosSims(i,22) = getCosSimilarity(adFCs(:,:,i), smad2FCs(:,:,i));
        cosSims(i,23) = getCosSimilarity(adFCs(:,:,i), smad7FCs(:,:,i));
        cosSims(i,24) = getCosSimilarity(adFCs(:,:,i), smad8FCs(:,:,i));
    end
    figure; boxplot(cosSims);
    
    % change Z score
%{
    cnDLWs = calcZScores(cnDLWs);
    adDLWs = calcZScores(adDLWs);
%}

    % normality test
%{
    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cnec', 'dlw');
%}
    % compalizon test (Wilcoxon, Mann?Whitney U test)
    [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcnFCs, roiNames, 'cn', 'smcn', 'fc');
    [adsmadFCsUt, adsmadFCsUtP, adsmadFCsUtP2] = calculateAlzWilcoxonTest(adFCs, smadFCs, roiNames, 'ad', 'smad', 'fc');
%    [cnsmrccnFCsUt, cnsmrccnFCsUtP, cnsmrccnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smrccnFCs, roiNames, 'cn', 'smrccn', 'fc');
%    [adsmrcadFCsUt, adsmrcadFCsUtP, adsmrcadFCsUtP2] = calculateAlzWilcoxonTest(adFCs, smrcadFCs, roiNames, 'ad', 'smrcad', 'fc');
    [cnsmcn2FCsUt, cnsmcn2FCsUtP, cnsmcn2FCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcn2FCs, roiNames, 'cn', 'smcn2', 'fc');
    [adsmad2FCsUt, adsmad2FCsUtP, adsmad2FCsUtP2] = calculateAlzWilcoxonTest(adFCs, smad2FCs, roiNames, 'ad', 'smad2', 'fc');
    [cnsmcn6FCsUt, cnsmcn6FCsUtP, cnsmcn6FCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcn6FCs, roiNames, 'cn', 'smcn6', 'fc');
    [cnsmcn7FCsUt, cnsmcn7FCsUtP, cnsmcn7FCsUtP2] = calculateAlzWilcoxonTest(cnFCs, smcn7FCs, roiNames, 'cn', 'smcn7', 'fc');
%    [cnsmcnFCsUt, cnsmcnFCsUtP, cnsmcnFCsUtP2] = calculateAlzWilcoxonTest(cnFCs, wtcnFCs, roiNames, 'cn', 'wtcn', 'fc');
    [cnsmcnDLsUt, cnsmcnDLsUtP, cnsmcnDLsUtP2] = calculateAlzWilcoxonTest(cnDLs, smcnDLs, roiNames, 'cn', 'smcn', 'dlcm');
    [adsmadDLsUt, adsmadDLsUtP, adsmadDLsUtP2] = calculateAlzWilcoxonTest(adDLs, smadDLs, roiNames, 'ad', 'smad', 'dlcm');
    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smcnDLWs, roiNames, 'cn', 'smcn', 'dlw');
    [adsmadDLWsUt, adsmadDLWsUtP, adsmadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smadDLWs, roiNames, 'ad', 'smad', 'dlw');
    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smcn7DLWs, roiNames, 'cn', 'smcn7', 'dlw');
%    [cnsmcnDLWsUt, cnsmcnDLWsUtP, cnsmcnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, wtcnDLWs, roiNames, 'cn', 'wtcn', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigCnDLWs, sigSmcnDLWs, roiNames, 'sigcn', 'sigsmcn', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigAdDLWs, sigSmadDLWs, roiNames, 'sigad', 'sigsmad', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigCnDLWs, sigSmcn7DLWs, roiNames, 'sigcn', 'sigsmcn7', 'dlw');
    [~, ~, ~] = calculateAlzWilcoxonTest(sigAdDLWs, sigSmad7DLWs, roiNames, 'sigad', 'sigsmad7', 'dlw');
%    [cnsmrccnDLWsUt, cnsmrccnDLWsUtP, cnsmrccnDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, smrccnDLWs, roiNames, 'cn', 'smrccn', 'dlw');
%    [adsmrcadDLWsUt, adsmrcadDLWsUtP, adsmrcadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, smrcadDLWs, roiNames, 'ad', 'smrcad', 'dlw');
end

% ==================================================================================================================
function [score, explained, pvals] = checkPCAHCvsAD(cnECs, adECs, nodeNum, algorithm)
    cnSbjNum = size(cnECs,3);
    adSbjNum = size(adECs,3);
    sbjECRoiVec = nan(cnSbjNum+adSbjNum,nodeNum*nodeNum-nodeNum);
    NE = eye(nodeNum); NE(NE==1) = nan;
    for i=1:cnSbjNum
        cnEC = cnECs(:,:,i) + NE;
        cnEC(isnan(cnEC)) = [];
        sbjECRoiVec(i,:) = cnEC(:);
    end
    for i=1:adSbjNum
        adEC = adECs(:,:,i) + NE;
        adEC(isnan(adEC)) = [];
        sbjECRoiVec(cnSbjNum+i,:) = adEC(:);
    end
    [coeff,score,latent,tsquared,explained,mu] = pca(sbjECRoiVec);
    % find well separated component by u test
    pvals = nan(size(score,2),1);
    for i=1:size(score,2)
        [pvals(i), h] = ranksum(score(1:cnSbjNum,i),score(cnSbjNum+1:end,i));
    end
    [B,I] = sort(pvals);
    figure; hold on;
    scatter3(score(1:cnSbjNum,I(1)),score(1:cnSbjNum,I(2)),score(1:cnSbjNum,I(3)),18,'r','filled');
    scatter3(score(cnSbjNum+1:end,I(1)),score(cnSbjNum+1:end,I(2)),score(cnSbjNum+1:end,I(3)),18,'b','filled');
    hold off; title(['PCA analysis of AD and HC regional vectors (' algorithm ')']); view(40,35)
    xlabel('component1'); ylabel('component2'); 
end

function [shiftDLWs, shiftSubDLWs, shiftSignals] = shiftAndExpandAmplitude(signals, DLWs, subDLWs, simSignals, simDLWs, simSubDLWs, group, type)
    nodeNum = size(simSignals{1},1);
    sigLen = size(simSignals{1},2);
    sbjNum = length(simSignals);
    R = nodeNum;
    nMax = 1;
    sbjMax = sbjNum;

    sfName = ['results/adsim2-shiftAmp' num2str(type) '-' group '-all.mat'];
    if exist(sfName, 'file')
        load(sfName);
        return;
    end

    if type == 2, nMax = 2; end
    if type == 3, amps = [2,4,6,8,10,12,14]; nMax = length(amps); end
    if type == 4, amps = [2,3,4,5,6,7,8]; nMax = length(amps); end
    if type == 5, nMax = 1; end

    shiftDLWs = simDLWs;
    shiftSubDLWs = simSubDLWs;
    shiftSignals = cell(sbjNum,1);
    cosSim = nan(sbjMax, nMax+1);

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 11;

    if NumProcessors > 1
        try
            disp('Destroing any existance matlab pool session');
            parpool('close');
        catch
            disp('No matlab pool session found');
        end
        parpool(NumProcessors);
    end

    % checking signal parallel shift effect for Zi, Zij and ECij'
    for k=1:sbjMax
        subEC = subDLWs(:,:,k);
        smSi = simSignals{k};
        smSubEC = simSubDLWs(:,:,k);
        sftSubEC = smSubEC;

        % calc & plot correlation of original Zi vs simulating Zi
        [corrZi, corrZij] = calcCorrelationZiZij(subEC, smSubEC, R);
        plotCorrelationZiZij([], subEC, [], smSubEC, R, ['sbj' num2str(k)], 'original', 'simulating');
        % plot original signals
%        figure; hold on; plot(smSi','Color',[0.8, 0.8, 0.8]); plot(smSi(1:R,:)');
%        hold off; title(['sbj' num2str(k) ' simulating signals']);

        outfName = ['results/adsim2-shiftAmp' num2str(type) '-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            sftSubECs = zeros(nodeNum,nodeNum+1,nMax);
            sftSignals = cell(1,nMax);
            corrZi2 = nan(1,nMax);
            corrZij2 = nan(R,nMax);
        end

        for n=1:nMax
            if type == 1
                % shift mean value of simulating signal
                if n >= 2, continue; end
                dZi = sftSubEC(:,1) - subEC(:,1); % smZi - Zi
                smSi = smSi(:,:) - dZi * 2 / 3;
                smSi(smSi>1.2) = 1.2;
                smSi(smSi<-0.2) = -0.2;
            elseif type == 2
                if n==1
                    amp = 2;
                    mvsi = movmean(smSi,65,2);
                    for i=1:R
                        smSi(i,:) = (smSi(i,:)-mvsi(i,:)) .* amp + mvsi(i,:);
                        smSi(i,:) = smSi(i,:) - nanmean(smSi(i,:));
                    end
                else
                    dZi = sftSubEC(:,1) - subEC(:,1); % smZi - Zi
                    smSi = smSi(:,:) - dZi * 2 / 3;
                    smSi(smSi>1.2) = 1.2;
                    smSi(smSi<-0.2) = -0.2;
                end
            elseif type == 3
                for i=1:nodeNum
                    wt=cwt(smSi(i,:));
                    wt(1:3,:) = wt(1:3,:) * amps(n);
                    trend = smoothdata(smSi(i,:),'movmean',32);
                    smSi(i,:) = icwt(wt,'SignalMean',trend);
                end
            elseif type == 4
                for i=1:nodeNum
                    wt=cwt(smSi(i,:));
                    wt(1:10,:) = wt(1:10,:) * amps(n);
                    trend = smoothdata(smSi(i,:),'movmean',32);
                    smSi(i,:) = icwt(wt,'SignalMean',trend);
                end
            elseif type == 5
                si = convert2SigmoidSignal(signals{k});
                for i=1:nodeNum
                    wtSi=cwt(si(i,:));
                    wtSm=cwt(smSi(i,:));
                    swtSi=nanstd(abs(wtSi),1,2);
                    swtSm=nanstd(abs(wtSm),1,2);
                    ml = swtSi ./ swtSm;
                    wtSm = wtSm .* ml;
                    trend = smoothdata(smSi(i,:),'movmean',32);
                    smSi(i,:) = icwt(wtSm,'SignalMean',trend);
%                    figure; cwt(smSi(i,:)'); figure; cwt(si(i,:)');
%                    figure; plot(smSi(i,:)'); figure; plot(si(i,:)');
                end
            end
            % amplitude expansion of simulating signal
            %{
            if n == 0
                alpha = 1;
                for i=1:R
                    sdsftZij = nanstd(sftSubEC(i,2:end),1);
                    sdZij = nanstd(subEC(i,2:end),1);
                    b = sdsftZij - 0.032;
                    amp = (sdZij - b) / 0.032 * alpha;
                    mvsi = movmean(smSi,5,2);
                    tsi = (smSi(i,:)-mvsi(i,:)) .* amp + mvsi(i,:);
                    tsi(tsi>1.2) = 1.2;
                    tsi(tsi<-0.2) = -0.2;
                    smSi(i,:) = tsi;
                end
%                smSi = convert2SigmoidSignal(smSi);
            end
            %}
            sftSignals{n} = smSi;

            % plot shifted simulating signals
%            figure; hold on; plot(smSi','Color',[0.8, 0.8, 0.8]); plot(smSi(1:R,:)');
%            hold off; title(['sbj' num2str(k) ' shifted simulating signals']);

            if isnan(corrZi2(n))
                % load original signal
                dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
                f = load(dlcmName);
                if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
                if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility

                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

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

    %            for i=1:R
                parfor i=1:R
                    netDLCM = initDlcmNetwork(smSi, f.exSignal, [], f.exControl); 

                    nodeTeach = smSi(i,2:end);
                    nodeInput = [smSi(:,1:end-1); f.exSignal(:,1:end-1)];
                    filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    disp(['training sbj' num2str(k) ' ' num2str(n) '-' num2str(i)]);
                    [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                    % predict DLCM network
                    sftSubEC(i,:) = predict(nodeNetwork, Si1);
                end
                sftSubECs(:,:,n) = sftSubEC;

                % calcate correlation
                [corrZi2(n), corrZij2(:,n)] = calcCorrelationZiZij(subEC, sftSubEC, R);

                save(outfName, 'corrZi', 'corrZi2', 'corrZij', 'corrZij2', 'sftSubECs', 'sftSignals');
            else
                sftSubEC = sftSubECs(:,:,n);
            end
            
            % plot correlation of original Zi vs shifted simulating Zi
%            plotCorrelationZiZij([], subEC, [], sftSubEC, R, ['sbj' num2str(k)], 'original', 'shifted sim');

            % cos similality
            sftEC = abs(repmat(sftSubEC(:,1),[1 nodeNum]) - sftSubEC(:,2:end));
            cosSim(k,n+1) = getCosSimilarity(DLWs(:,:,k), sftEC);
        end

        % find most correlated Zi & Zij signals
        cosSim(k,1) = getCosSimilarity(DLWs(:,:,k), simDLWs(:,:,k));
        [m, idx] = max(cosSim(k,:));
        if idx==1
            sftSignal = simSignals{k};
            sftSubEC = simSubDLWs(:,:,k);
        else
            sftSignal = sftSignals{idx-1};
            sftSubEC = sftSubECs(:,:,idx-1);
        end

        % plot correlation of original Zi vs shifted simulating Zi
%        plotCorrelationZiZij([], subEC, [], sftSubEC, R, ['sbj' num2str(k)], 'original', 'shifted sim');

        % set output matrix
        shiftDLWs(:,:,k) = abs(repmat(sftSubEC(:,1),[1 nodeNum]) - sftSubEC(:,2:end));
        shiftSubDLWs(:,:,k) = sftSubEC;
        shiftSignals{k} = sftSignal;
    end
    save(sfName, 'shiftDLWs', 'shiftSubDLWs', 'shiftSignals', 'cosSim');

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
end

function checkRelationSubDLWandWeights(signals, subDLWs, smSubDLWs, group, type)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    R = 1; %nodeNum;
    sbjMax = 2;

    % checking signal parallel shift effect for Zi, Zij and ECij'
    for k=1:sbjMax
        Zi = subDLWs(:,1,k);
        Zij = subDLWs(:,2:end,k);
        smZi = smSubDLWs(:,1,k);
        smZij = smSubDLWs(:,2:end,k);
        ECd = Zij - repmat(Zi,[1,nodeNum]);
        smECd = smZij - repmat(smZi,[1,nodeNum]);
        
        outfName = ['results/adsim2-checkRelation7-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            
            % training options for DLCM network
            options = trainingOptions('adam', 'InitialLearnRate', 0.0001, 'ExecutionEnvironment','cpu', 'MaxEpochs',1, 'Verbose',false);
                
            % training loop
            for i=1:R 
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                % change weight value
                tmpL = f.netDLCM.nodeNetwork{i}.Layers;
                for j=1:nodeNum
                    if i==j, continue; end
                    if type == 1
                        dx = ECd(i,j) - smECd(i,j);
                        tmpL(2).Weights(:,j) = tmpL(2).Weights(:,j) - dx * 0.7;
                    elseif type == 2
                        dx = smECd(i,j) / ECd(i,j);
                        tmpL(2).Weights(:,j) = tmpL(2).Weights(:,j) / dx;
                    end
                end
                Layers = makeDlcmLayers(tmpL);

                nodeTeach = f.si(i,2:end);
                nodeInput = [f.si(:,1:end-1); f.exSignal(:,1:end-1)];
                filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                disp(['training ' num2str(k) '-' num2str(i)]);
                [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, Layers, options);
                % predict DLCM network
                smSubEC2 = predict(nodeNetwork, Si1);
                smZi2 = smSubEC2(1);
                smZij2 = smSubEC2(2:end);
                smECd2 = smZij2 - repmat(smZi2,[1,nodeNum]);
%%{
                for j=3:3
                    % plot Zi-Zij scat
                    figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                    scatter(ECd(i,:),smECd(i,:),3,[0.7 0.7 0.7]);
                    scatter(ECd(i,:),smECd2(:),3,[0.4 0.3 0.3]);
                    scatter(ECd(i,j),smECd2(j),3,[0.7 0.2 0.2]);
                    hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : original vs shifted sim']);

                    % plot Zi-Zij scat
                    figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                    scatter(smECd(i,:),smECd2(:),3,[0.3 0.3 0.3]);
                    scatter(smECd(i,j),smECd2(j),3,[0.7 0.2 0.2]);
                    hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : simulated vs shifted sim']);
                end
%%}
            end
        end
    end
end

function [weDLWs, weSubDLWs, weSignals, weDLs] = checkRelationSubDLWandWeights2(signals, DLWs, subDLWs, smSignals, smDLWs, smSubDLWs, smDLs, group, type)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    sbjMax = sbjNum;

    typename = '';
    if type==2, typename='t2'; end
    
    sfName = ['results/adsim2-checkRelation8' typename '-' group '-all.mat'];
    if exist(sfName, 'file')
        load(sfName);
        return;
    end
    
    weDLWs = DLWs;
    weSubDLWs = subDLWs;
    weSignals = cell(sbjNum,1);
    weDLs = nan(nodeNum,nodeNum,sbjNum);

    if type == 1
        wes = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
    elseif type == 2
        wes = [0.1,0.3,1];
    end
    weLen = length(wes);

    % calc correlation
    ECdCr = nan(sbjMax, nodeNum);
    smECd2Cr = nan(sbjMax, weLen+1, nodeNum);
    cosSim = nan(sbjMax, weLen+1,1);

    % checking signal parallel shift effect for Zi, Zij and ECij'
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        smEC = smDLWs(:,:,k);
        Zi = subDLWs(:,1,k);
        Zij = subDLWs(:,2:end,k);
        smZi = smSubDLWs(:,1,k);
        smZij = smSubDLWs(:,2:end,k);
        ECd = Zij - repmat(Zi,[1,nodeNum]);
        smECd = smZij - repmat(smZi,[1,nodeNum]);
        
        outfName = ['results/adsim2-checkRelation8' typename '-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            tmpfName = ['results/adsim2-checkRelation8' typename '-' group '-' num2str(k) '_tmp.mat'];
            if exist(tmpfName, 'file')
                load(tmpfName);
                westart = length(weSi) + 1;
            else
                EC2s = cell(weLen,1);
                subEC2s = cell(weLen,1);
                weSi = {};
                weSi2 = {};
                weNet = {};
                weNet2 = {};
                westart = 1;
            end

            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility

            % training options for weight expansion training
            weOptions = trainingOptions('adam', 'InitialLearnRate', 0.0001, 'ExecutionEnvironment','cpu', 'MaxEpochs',1, 'Verbose',false);
            % training options for DLCM network
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

            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 11;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end

            netDLCM = f.netDLCM;
            exSignal = f.exSignal;
            exControl = f.exControl;
            [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{k});
            weSi{end+1} = si;

            for a=westart:weLen
                weRate = wes(a);
                nodeNetwork = cell(nodeNum,1);

                % weight expansion training loop
%                for i=1:nodeNum
                parfor i=1:nodeNum
                    % change weight value
                    tmpL = f.netDLCM.nodeNetwork{i}.Layers;
                    for j=1:nodeNum
                        if i==j, continue; end
                        if type == 1
                            dx = ECd(i,j) - smECd(i,j);
                            tmpL(2).Weights(:,j) = tmpL(2).Weights(:,j) - dx * weRate;
                        elseif type == 2
                            dx = smECd(i,j) / ECd(i,j);
                            tmpL(2).Weights(:,j) = tmpL(2).Weights(:,j) / dx * weRate;
                        end
                    
                    end
                    layers = makeDlcmLayers(tmpL);

                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); exSignal(:,1:end-1)];
                    filter = repmat(exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    disp(['sbj' num2str(k) ' weight expansion training ' num2str(i)]);
                    [nodeNetwork{i}, ~] = trainNetwork(nodeInput, nodeTeach, layers, weOptions);
                end
                netDLCM.nodeNetwork = nodeNetwork;
                weNet{end+1} = netDLCM;

                % simulate signal from first frame
                [si2, time] = simulateDlcmNetwork(si, exSignal, [], exControl, netDLCM);
                weSi2{end+1} = si2;

                % train DLCM network with simulated signal of amplitude expanded DLCM network
                netDLCM2 = initDlcmNetwork(si2, exSignal, [], exControl); 

                disp(['sbj' num2str(k) ' training 2nd']);
                netDLCM2 = trainDlcmNetwork(si2, exSignal, [], exControl, netDLCM2, options);
                weNet2{end+1} = netDLCM2;

                % calculate DLCM-EC
                [EC2s{a}, subEC2s{a}] = calcDlcmEC(netDLCM2, [], exControl);

                save(tmpfName, 'EC2s', 'subEC2s', 'exSignal', 'exControl', 'weSi', 'weSi2', 'weNet', 'weNet2');
            end

            % calc correlation between ECd and smECd2
            for a=1:weLen
                weRate = wes(a);
                smZi2 = subEC2s{a}(:,1);
                smZij2 = subEC2s{a}(:,2:end);
                smECd2 = smZij2 - repmat(smZi2,[1,nodeNum]);

                for i=1:nodeNum
                    ECdCr(k,i) = ecCorr(ECd(i,:),smECd(i,:), i);
                    smECd2Cr(k,a+1,i) = ecCorr(ECd(i,:),smECd2(i,:), i);

                    % plot Zi-Zij scat
%                    figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
%                    scatter(ECd(i,:),smECd(i,:),3,[0.8 0.8 0.8]);
%                    scatter(ECd(i,:),smECd2(:),3,[0.5 0.3 0.3]);
%                    hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : original vs shifted sim']);
                end

                % plot Zi-Zij all node
%{
                figure; hold on; plot([-0.6 0.6], [-0.6 0.6],':','Color',[0.5 0.5 0.5]);
                for i=1:nodeNum
                    plotTwoSignalsCorrelation(ECd(i,:),smECd(i,:), [0.4+0.2*ceil(i/100) 0.1*mod(i,5) 0.1*ceil(mod(i,50)/10)], 'x', 5);
                end
                hold off; title(['sbj' num2str(k) ' rate=' num2str(weRate) ' Zij-Zi corr : original vs simulated']);
                figure; hold on; plot([-0.6 0.6], [-0.6 0.6],':','Color',[0.5 0.5 0.5]);
                for i=1:nodeNum
                    plotTwoSignalsCorrelation(ECd(i,:),smECd2(i,:), [0.4+0.2*ceil(i/100) 0.1*mod(i,5) 0.1*ceil(mod(i,50)/10)], 'x', 5);
                end
                hold off; title(['sbj' num2str(k) ' rate=' num2str(weRate) ' Zij-Zi corr : original vs shifted sim']);
                % plot Zi & Zij all node
                plotCorrelationZiZij([], subDLWs(:,:,k), [], subEC2s{a}, nodeNum, ['sbj' num2str(k) ' rate=' num2str(weRate)], 'original', 'shifted sim');
%}
            end
            save(outfName, 'EC2s', 'subEC2s', 'ECdCr', 'smECd2Cr', 'exSignal', 'exControl', 'weSi', 'weSi2', 'weNet', 'weNet2');
            
            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
        
        % calc cos similarity
        cosSim(k,1) = getCosSimilarity(EC, smEC);
        for a=1:weLen
            cosSim(k,a+1) = getCosSimilarity(EC, EC2s{a});
        end
        % find most similar signal
        [m,idx] = max(cosSim(k,:));
        if idx==1
            weSignals{k} = smSignals{k};
            weDLWs(:,:,k) = smDLWs(:,:,k);
            weSubDLWs(:,:,k) = smSubDLWs(:,:,k);
            weDLs(:,:,k) = smDLs(:,:,k);
        else
            weSignals{k} = weSi2{idx-1};
            weDLWs(:,:,k) = EC2s{idx-1};
            weSubDLWs(:,:,k) = subEC2s{idx-1};
            weDLs(:,:,k) = calcDlcmGCI(weSi2{idx-1}, exSignal, [], exControl, weNet2{idx-1});
        end
    end
    save(sfName, 'weDLWs', 'weSubDLWs', 'weSignals', 'weDLs', 'ECdCr', 'smECd2Cr', 'cosSim');
end

function cr = ecCorr(X, Y, i)
    X(i) = nan;
    Y(i) = nan;
    A = X(~isnan(X));
    B = Y(~isnan(Y));
    cr = corr2(A(:), B(:));
end

function layers = makeDlcmLayers(oldLayers)
    % init first fully connected layer
    inNum = size(oldLayers(2).Weights,2);
    hiddenNums(1) = size(oldLayers(2).Weights,1);
    hiddenNums(2) = size(oldLayers(4).Weights,1);
    
    inLayers = [
        % input layer
        sequenceInputLayer(inNum);
        % Add a fully connected layer
        fullyConnectedLayer(hiddenNums(1), 'Weights', oldLayers(2).Weights, 'Bias', oldLayers(2).Bias);
        % Add an ReLU non-linearity.
        reluLayer();
        ];
    hdLayers = [
        % Add a fully connected layer
        fullyConnectedLayer(hiddenNums(2), 'Weights', oldLayers(4).Weights, 'Bias', oldLayers(4).Bias)
        % Add an ReLU non-linearity.
        reluLayer();
    ];
    layers = [
        inLayers;
        hdLayers;
        % Add a fully connected layer
        fullyConnectedLayer(1, 'Weights', oldLayers(6).Weights, 'Bias', oldLayers(6).Bias);
        % reggression for learning
        regressionLayer();
    ];
end

function checkRelationSubDLWandSignals4(signals, DLWs, subDLWs, smSignals, smDLWs, smSubDLWs, group)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    R = 1; %nodeNum;
    nMax = 1;
    sbjMax = 1;

    % checking signal parallel shift effect for Zi, Zij and ECij'
    for k=1:sbjMax
        Zi = subDLWs(:,1,k);
        Zij = subDLWs(:,2:end,k);
        smZi = smSubDLWs(:,1,k);
        smZij = smSubDLWs(:,2:end,k);
        ECd = Zij - repmat(Zi,[1,nodeNum]);
        smECd = smZij - repmat(smZi,[1,nodeNum]);
        
        outfName = ['results/adsim2-checkRelation7-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            for i=1:R
                % plot Zi-Zij scat
                figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                scatter(ECd(i,:),smECd(i,:),3,[0.3 0.3 0.3]);
                scatter(ECd(i,2),smECd(i,2),3,[0.7 0.2 0.2]);
                hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : original vs simulated']);
            end

            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            siOrg = smSignals{k};
            
            % training options for DLCM network
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
                
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 1;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end

            % training loop
            for i=1:R 
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                si = siOrg;
                for j=1:nodeNum
                    if i==j, continue; end
                    % shift target signal
                    % -- parallel does not affect
                    dx = ECd(i,j) - smECd(i,j);
                    si(j,:) = si(j,:) + dx;
                    % -- amplitude expansion affect badly
%                    amp = 2;
%                    mvsi = movmean(si,65,2);
%                    si(j,:) = (si(j,:)-mvsi(j,:)) .* amp + mvsi(j,:);
%                    si(j,:) = si(j,:) + (nanmean(siOrg(j,:)) - nanmean(si(j,:)));

                    % plot original & shifted signal
%                    figure; hold on; plot(siOrg','Color',[0.7 0.7 0.7]); plot(siOrg(j,:)','Color',[0.2 0.2 0.7]); plot(si(j,:)','Color',[0.7 0.2 0.2]); hold off;
                end
                figure; hold on; plot(siOrg','Color',[0.7 0.7 0.7]); plot(si','Color',[0.8 0.4 0.4]); hold off;

                % DLCM network training
                netDLCM = initDlcmNetwork(si, f.exSignal, [], f.exControl); 

                nodeTeach = si(i,2:end);
                nodeInput = [si(:,1:end-1); f.exSignal(:,1:end-1)];
                filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                disp(['training ' num2str(k) '-' num2str(i)]);
                [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                % predict DLCM network
                smSubEC2 = predict(nodeNetwork, Si1);
                smZi2 = smSubEC2(1);
                smZij2 = smSubEC2(2:end);
                smECd2 = smZij2 - repmat(smZi2,[1,nodeNum]);
%{
                for j=1:nodeNum
                    % plot Zi-Zij scat
                    figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                    scatter(ECd(i,:),smECd2(:),3,[0.3 0.3 0.3]);
                    scatter(ECd(i,j),smECd2(j),3,[0.7 0.2 0.2]);
                    hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : original vs shifted sim']);

                    % plot Zi-Zij scat
                    figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                    scatter(smECd(i,:),smECd2(:),3,[0.3 0.3 0.3]);
                    scatter(smECd(i,j),smECd2(j),3,[0.7 0.2 0.2]);
                    hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : simulated vs shifted sim']);
                end
%}
                % plot Zi-Zij scat
                figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                scatter(ECd(i,:),smECd(i,:),3,[0.7 0.7 0.7]);
                scatter(ECd(i,:),smECd2(:),3,[0.7 0.2 0.2]);
                hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : original vs shifted sim']);                
                % plot Zi-Zij scat
                figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                scatter(smECd(i,:),smECd2(:),3,[0.3 0.3 0.3]);
                hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : simulated vs shifted sim']);
            end

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
    end
end

function [wtDLWs, wtSubDLWs, wtSignals, wtDLs] = checkWaveletTransformEffect(signals, DLWs, subDLWs, group, type)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    sbjMax = sbjNum;

    typename = '';
    if type==2, typename='t2'; end
    if type==3, typename='t3'; end
    if type==4, typename='t4'; end

    sfName = ['results/adsim2-checkWavelet' typename '-' group '-all.mat'];
    if exist(sfName, 'file')
        load(sfName);
        return;
    end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 11;

    if NumProcessors > 1
        try
            disp('Destroing any existance matlab pool session');
            parpool('close');
        catch
            disp('No matlab pool session found');
        end
        parpool(NumProcessors);
    end

    wtDLWs = DLWs;
    wtSubDLWs = subDLWs;
    wtSignals = cell(sbjNum,1);
    wtDLs = nan(nodeNum,nodeNum,sbjNum);

    % training options for DLCM network
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
        
    % checking signal amplitude change effect for Zi, Zij and ECij'
    for k=1:sbjMax
        dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
        f = load(dlcmName);
        if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
        if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
        exSignal = f.exSignal;
        exControl = f.exControl;

        % wavelet transform and back
        siOrg = signals{k};
        if type==1
            for i=1:nodeNum
                wt=cwt(siOrg(i,:));
                trend = smoothdata(siOrg(i,:),'movmean',32);
                siOrg(i,:) = icwt(wt,'SignalMean',trend);
            end
        elseif type==2
            for i=1:nodeNum
                [S,F,T] = stft(siOrg(i,:),'Window',hamming(24,'periodic'),'OverlapLength',16,'FFTLength',64);
                siOrg(i,:) = istft(S,'Window',hamming(24,'periodic'),'OverlapLength',16,'FFTLength',64);
            end
        elseif type==3
            for i=1:nodeNum
                wt=cwt(siOrg(i,:));
                wt(1:15,:) = 0;
                trend = smoothdata(siOrg(i,:),'movmean',32);
                siOrg(i,:) = icwt(wt,'SignalMean',trend);
            end
        elseif type==4
            for i=1:nodeNum
                wt=cwt(siOrg(i,:));
                siOrg(i,:) = icwt(wt);
            end
        end
        [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(siOrg);
        wtSignals{k} = si;

        % plot original and modfied
        figure; plot(si'); figure; si2 = convert2SigmoidSignal(signals{k}); plot(si2'); title(['subject' num2str(k)]);

        % train DLCM network with amplitude expanded signal
        netDLCM = initDlcmNetwork(si, exSignal, [], exControl); 

        disp(['sbj' num2str(k) ' training']);
        netDLCM = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);

        % calculate DLCM-EC
        [wtDLWs(:,:,k), wtSubDLWs(:,:,k)] = calcDlcmEC(netDLCM, [], exControl);
        wtDLs(:,:,k) = calcDlcmGCI(si, exSignal, [], exControl, netDLCM);
    end
    save(sfName, 'wtDLWs', 'wtSubDLWs', 'wtSignals', 'wtDLs');

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
end

function [ampDLWs, ampSubDLWs, ampSignals, ampDLs] = checkRelationSubDLWandSignals3(signals, DLWs, subDLWs, smSignals, smDLWs, smSubDLWs, smDLs, group, type)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    sbjMax = sbjNum;

    typename = '';
    if type==2, typename='t2'; end
    if type==3, typename='t3'; end
    if type==4, typename='t4'; end
    if type==5, typename='t5'; end
    if type==6, typename='t6'; end
    
    sfName = ['results/adsim2-checkRelation6' typename '-' group '-all.mat'];
    if exist(sfName, 'file')
        load(sfName);
        return;
    end
    
    ampDLWs = DLWs;
    ampSubDLWs = subDLWs;
    ampSignals = cell(sbjNum,1);
    ampDLs = nan(nodeNum,nodeNum,sbjNum);

    if type == 1
        amps = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
    elseif type == 2
        amps = [0.5, 1, 1.5, 2, 2.5];
    elseif type == 3
        amps = [0.5, 1, 1.5];
    elseif type == 4
        amps = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
    elseif type == 5
        amps = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
    elseif type == 6
        amps = [2, 4, 6, 8, 10, 12, 14];
    end
    ampsLen = length(amps);

    % calc correlation
    ZiCr = nan(sbjMax, ampsLen+1, 1);
    ZijCr = nan(sbjMax, ampsLen+1, nodeNum);
    cosSim = nan(sbjMax, ampsLen+1,1);

    % checking signal amplitude change effect for Zi, Zij and ECij'
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);
        smEC = smDLWs(:,:,k);
        smSubEC = smSubDLWs(:,:,k);

        Zi = subDLWs(:,1,k);
        Zij = subDLWs(:,2:end,k);
        smZi = smSubDLWs(:,1,k);
        smZij = smSubDLWs(:,2:end,k);
        ECd = Zij - repmat(Zi,[1,nodeNum]);
        smECd = smZij - repmat(smZi,[1,nodeNum]);

        outfName = ['results/adsim2-checkRelation6' typename '-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            tmpfName = ['results/adsim2-checkRelation6' typename '-' group '-' num2str(k) '_tmp.mat'];
            if exist(tmpfName, 'file')
                load(tmpfName);
                astart = length(ampSi) + 1;
            else
                EC2s = cell(ampsLen,1);
                subEC2s = cell(ampsLen,1);
                ampSi = {};
                ampSi2 = {};
                ampNet = {};
                ampNet2 = {};
                astart = 1;
            end
            
            % plot Zi-Zij scat
%{
            for j=1:nodeNum
                figure; hold on; plot([-0.2 0.2], [-0.2 0.2],':','Color',[0.5 0.5 0.5]);
                scatter(ECd(:),smECd(:),3,[0.7 0.7 0.7]);
                scatter(ECd(j,:),smECd(j,:),3,[0.2 0.2 0.7]);
                scatter(ECd(:,j),smECd(:,j),3,[0.7 0.2 0.2]);
                hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij-Zi corr : original vs shifted sim']);
            end
%}
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 11;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            exSignal = f.exSignal;
            exControl = f.exControl;

            % training options for DLCM network
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

            for a=astart:ampsLen
                amp = amps(a);
                siOrg = signals{k};
                if type == 1
                    mvsi = movmean(siOrg,65,2);
                    for i=1:nodeNum % target node
                        siOrg(i,:) = (siOrg(i,:)-mvsi(i,:)) .* amp + mvsi(i,:);
                        siOrg(i,:) = siOrg(i,:) - nanmean(siOrg(i,:));
                    end
                elseif type == 2
                    orgStd = nanstd(ECd,1);
                    smStd = nanstd(smECd,1);
                    dAmp = orgStd ./ smStd;
                    mvsi = movmean(siOrg,65,2);
                    for i=1:nodeNum % target node
                        siOrg(i,:) = (siOrg(i,:)-mvsi(i,:)) * dAmp(i) * amp + mvsi(i,:);
                        siOrg(i,:) = siOrg(i,:) - nanmean(siOrg(i,:));
                    end
                elseif type == 3
                    orgStd = nanstd(ECd,1,2);
                    smStd = nanstd(smECd,1,2);
                    dAmp = orgStd ./ smStd;
                    mvsi = movmean(siOrg,65,2);
                    for i=1:nodeNum % target node
                        siOrg(i,:) = (siOrg(i,:)-mvsi(i,:)) * dAmp(i) * amp + mvsi(i,:);
                        siOrg(i,:) = siOrg(i,:) - nanmean(siOrg(i,:));
                    end
                elseif type == 4
                    %{
                    xDFT = fft(siOrg');
                    figure; plot(siOrg'); xlabel('Seconds'); ylabel('Amplitude');
                    sz = size(xDFT,1)/2 + 1;
                    figure; plot(abs(xDFT(1:sz,:)));
                    set(gca,'xtick',[4:4:64]);
                    xlabel('Hz'); ylabel('Magnitude');
                    %}
                    for i=1:nodeNum
                        wt=cwt(siOrg(i,:));
                        wt(1:15,:)=wt(1:15,:) * amp;
                        trend = smoothdata(siOrg(i,:),'movmean',32);
                        siOrg(i,:) = icwt(wt,'SignalMean',trend);
                    end
                elseif type == 5
                    for i=1:nodeNum
                        wt=cwt(siOrg(i,:));
                        wt(31:end,:)=wt(31:end,:) * amp;
                        trend = smoothdata(siOrg(i,:),'movmean',32);
                        siOrg(i,:) = icwt(wt,'SignalMean',trend);
                    end
                elseif type == 6
                    for i=1:nodeNum
                        wt=cwt(siOrg(i,:));
                        wt(1:3,:)=wt(1:3,:) * amp;
                        trend = smoothdata(siOrg(i,:),'movmean',32);
                        siOrg(i,:) = icwt(wt,'SignalMean',trend);
                    end
                end
                [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(siOrg);
                ampSi{end+1} = si;
                
                % plot original and modfied
%                figure; plot(si'); figure; si2 = convert2SigmoidSignal(signals{k}); plot(si2');

                % train DLCM network with amplitude expanded signal
                netDLCM = initDlcmNetwork(si, exSignal, [], exControl); 

                disp(['sbj' num2str(k) ' training 1st amp=' num2str(amp)]);
                netDLCM = trainDlcmNetwork(si, exSignal, [], exControl, netDLCM, options);
                ampNet{end+1} = netDLCM;

                % simulate signal from first frame
                [si2, time] = simulateDlcmNetwork(si, exSignal, [], exControl, netDLCM);
                ampSi2{end+1} = si2;

                % train DLCM network with simulated signal of amplitude expanded DLCM network
                netDLCM2 = initDlcmNetwork(si2, exSignal, [], exControl); 

                disp(['sbj' num2str(k) ' training 2nd amp=' num2str(amp)]);
                netDLCM2 = trainDlcmNetwork(si2, exSignal, [], exControl, netDLCM2, options);
                ampNet2{end+1} = netDLCM2;

                % calculate DLCM-EC
                [EC2s{a}, subEC2s{a}] = calcDlcmEC(netDLCM2, [], exControl);

                save(tmpfName, 'EC2s', 'subEC2s', 'exSignal', 'exControl', 'ampSi', 'ampSi2', 'ampNet', 'ampNet2');
            end
            save(outfName, 'EC2s', 'subEC2s', 'exSignal', 'exControl', 'ampSi', 'ampSi2', 'ampNet', 'ampNet2');
            
            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end

        % plot & calc correlation of original vs simulating
        plotCorrelationZiZij(EC, subEC, smEC, smSubEC, nodeNum, ['sbj' num2str(k)], 'original', 'simulated');
        [ZiCr(k,1), ZijCr(k,1,:)] = calcCorrelationZiZij(subEC, smSubEC, nodeNum);

        % plot & calc correlation of original vs shifted sim
        for a=1:ampsLen
%            plotCorrelationZiZij(EC, subEC, EC2s{a}, subEC2s{a}, nodeNum, ['sbj' num2str(k) ' amp=' num2str(amps(a))], 'original', 'shifted sim');
            [ZiCr(k,a+1), ZijCr(k,a+1,:)] = calcCorrelationZiZij(subEC, subEC2s{a}, nodeNum);
        end

        % calc cos similarity
        cosSim(k,1) = getCosSimilarity(EC, smEC);
        for a=1:ampsLen
            cosSim(k,a+1) = getCosSimilarity(EC, EC2s{a});
        end
        % find most similar signal
        [m,idx] = max(cosSim(k,:));
        if idx==1
            ampSignals{k} = smSignals{k};
            ampDLWs(:,:,k) = smDLWs(:,:,k);
            ampSubDLWs(:,:,k) = smSubDLWs(:,:,k);
            ampDLs(:,:,k) = smDLs(:,:,k);
        else
            ampSignals{k} = ampSi2{idx-1};
            ampDLWs(:,:,k) = EC2s{idx-1};
            ampSubDLWs(:,:,k) = subEC2s{idx-1};
            ampDLs(:,:,k) = calcDlcmGCI(ampSi2{idx-1}, exSignal, [], exControl, ampNet2{idx-1});
        end
    end
    save(sfName, 'ampDLWs', 'ampSubDLWs', 'ampSignals', 'ampDLs', 'ZiCr', 'ZijCr', 'cosSim');
end

function [ampDLWs, ampSubDLWs, ampSignals, ampDLs] = checkRelationSubDLWandSignals3b(signals, DLWs, subDLWs, smSignals, smDLWs, smSubDLWs, smDLs, group, orgGroup)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    sbjMax = sbjNum;

    sfName = ['results/adsim2-checkRelation6b-' group '-all.mat'];
    if exist(sfName, 'file')
        load(sfName);
        return;
    end
    
    ampDLWs = DLWs;
    ampSubDLWs = subDLWs;
    ampSignals = cell(sbjNum,1);
    ampDLs = nan(nodeNum,nodeNum,sbjNum);

    amps = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5];
    ampsLen = length(amps);

    % calc correlation
    ZiCr = nan(sbjMax, ampsLen+1, 1);
    ZijCr = nan(sbjMax, ampsLen+1, nodeNum);
    cosSim = nan(sbjMax, ampsLen+1,1);

    % checking signal amplitude change effect for Zi, Zij and ECij'
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);
        smEC = smDLWs(:,:,k);
        smSubEC = smSubDLWs(:,:,k);
        
        outfName = ['results/adsim2-checkRelation6b-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            EC2s = cell(ampsLen,1);
            subEC2s = cell(ampsLen,1);
            ampSi2 = {};
            ampNet2 = {};
            
            ampfName = ['results/adsim2-checkRelation6-' orgGroup '-' num2str(k) '.mat'];
            af=load(ampfName);        

            dlcmName = ['results/ad-dlcm-' orgGroup '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            exSignal = f.exSignal;
            exControl = f.exControl;

            % and change for loop to parfor loop
            NumProcessors = 20;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end

            % training options for DLCM network
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

            for a=1:ampsLen
                % simulate signal from first frame
                [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{k});
                [si2, time] = simulateDlcmNetwork(si, exSignal, [], exControl, af.ampNet{a});
                ampSi2{end+1} = si2;

                % train DLCM network with simulated signal of amplitude expanded DLCM network
                netDLCM2 = initDlcmNetwork(si2, exSignal, [], exControl); 

                disp(['sbj' num2str(k) ' training 2nd amp=' num2str(amps(a))]);
                netDLCM2 = trainDlcmNetwork(si2, exSignal, [], exControl, netDLCM2, options);
                ampNet2{end+1} = netDLCM2;

                % calculate DLCM-EC
                [EC2s{a}, subEC2s{a}] = calcDlcmEC(netDLCM2, [], exControl);
            end

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
            save(outfName, 'EC2s', 'subEC2s', 'exSignal', 'exControl', 'ampSi2', 'ampNet2');
        end
        
        % plot & calc correlation of original vs simulating
        plotCorrelationZiZij(EC, subEC, smEC, smSubEC, nodeNum, ['sbj' num2str(k)], 'original', 'simulated');
        [ZiCr(k,1), ZijCr(k,1,:)] = calcCorrelationZiZij(subEC, smSubEC, nodeNum);

        % calc cos similarity & corrleation
        cosSim(k,1) = getCosSimilarity(EC, smEC);
        for a=1:ampsLen
            cosSim(k,a+1) = getCosSimilarity(EC, EC2s{a});

            % plot & calc correlation of original vs simulating
%            plotCorrelationZiZij(EC, subEC, EC2s{a}, subEC2s{a}, nodeNum, ['sbj' num2str(k) ' amp=' num2str(amps(a))], 'original', 'shifted sim');
            [ZiCr(k,a+1), ZijCr(k,a+1,:)] = calcCorrelationZiZij(subEC, subEC2s{a}, nodeNum);
        end
        % find most similar signal
        [m,idx] = max(cosSim(k,:));
        if idx==1
            ampSignals{k} = smSignals{k};
            ampDLWs(:,:,k) = smDLWs(:,:,k);
            ampSubDLWs(:,:,k) = smSubDLWs(:,:,k);
            ampDLs(:,:,k) = smDLs(:,:,k);
        else
            ampSignals{k} = ampSi2{idx-1};
            ampDLWs(:,:,k) = EC2s{idx-1};
            ampSubDLWs(:,:,k) = subEC2s{idx-1};
            ampDLs(:,:,k) = calcDlcmGCI(ampSi2{idx-1}, exSignal, [], exControl, ampNet2{idx-1});
        end
    end
    save(sfName, 'ampDLWs', 'ampSubDLWs', 'ampSignals', 'ampDLs', 'ZiCr', 'ZijCr', 'cosSim');
end

function plotCorrelationZiZij(EC, subEC, smEC, smSubEC, nodeNum, prefix, orig, sim)
    if ~isempty(subEC)
        figure; hold on; plot([0.5 1.2], [0.5 1.2],':','Color',[0.5 0.5 0.5]);
        for i=1:nodeNum
            X = subEC(i,2:end);
            Y = smSubEC(i,2:end); Y(i) = nan;
            plotTwoSignalsCorrelation(X, Y, [0.3+0.07*mod(i,10) 0.3+0.07*ceil(mod(i,100)/10) 0.6+0.2*ceil(i/100)]);
        end
        for i=1:nodeNum
            plotTwoSignalsCorrelation(subEC(i,1), smSubEC(i,1), [0.4+0.2*ceil(i/100) 0.1*mod(i,5) 0.1*ceil(mod(i,50)/10)], 'd', 8);
        end
        hold off; title([prefix ' Zij corr: ' orig ' vs ' sim]);
    end
    if ~isempty(EC)
        figure; hold on; plot([0 0.5], [0 0.5],':','Color',[0.5 0.5 0.5]);
        for i=1:nodeNum
            plotTwoSignalsCorrelation(EC(i,:), smEC(i,:), [0.4+0.2*ceil(i/100) 0.1*mod(i,5) 0.1*ceil(mod(i,50)/10)], 'x', 5);
        end
        hold off; title([prefix ' Zi corr: ' orig ' vs ' sim]);
    end
end

function [ZiCorr, ZijCorr] = calcCorrelationZiZij(subEC, smSubEC, nodeNum)        
    ZiCorr = corr2(subEC(:,1), smSubEC(:,1));
    ZijCorr = nan(nodeNum,1);
    for i=1:nodeNum
        X = subEC(i,2:end);
        Y = smSubEC(i,2:end); Y(i) = nan;
        ZijCorr(i) = corr2(X(~isnan(X)), Y(~isnan(Y)));
    end
end

function checkRelationSubDLWandSignals2(signals, DLWs, subDLWs, smDLWs, smSubDLWs, group)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    R = 4;
    nMax = 1;
    sbjMax = 1; %4;

    % checking signal amplitude change effect for Zi, Zij and ECij'
%%{
%    amps = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5, 2, 3, 5, 8];
    amps = [2, 3];
    ampsLen = length(amps);
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);
        smEC = smDLWs(:,:,k);
        smSubEC = smSubDLWs(:,:,k);
        
        outfName = ['results/adsim2-checkRelation4-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            tmpfName = ['results/adsim2-checkRelation4-' group '-' num2str(k) '_tmp.mat'];
            if exist(tmpfName, 'file')
                load(tmpfName);
                if size(X,1) < R
                    X(R, ampsLen, nMax) = 0;
                end
                idx = find(X(:, 1, 1)==0);
                rstart =idx(1);
            else
                Zi2 = zeros(R, ampsLen, nMax);
                X = zeros(R, ampsLen, nMax);
                Zij2 = zeros(R, ampsLen, nMax, nodeNum);
                subEC2s = cell(R,ampsLen,nMax);
                rstart = 1;
            end
            
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 11;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility

            % training options for DLCM network
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

            for i=rstart:R % target node
                for a=1:ampsLen
                    amp = amps(a);
                    siOrg = signals{k};
%%{
                    mvsi = movmean(siOrg,65,2);
                    siOrg(i,:) = (siOrg(i,:)-mvsi(i,:)) .* amp + mvsi(i,:);
                    siOrg(i,:) = siOrg(i,:) - nanmean(siOrg(i,:));
                    figure; hold on; plot(signals{k}(i,:)'); plot(mvsi(i,:)'); plot(siOrg(i,:)'); hold off; title(['node' num2str(i) ' amp=' num2str(amp)]);
%%}
%{
                    m = nanmean(siOrg(i,:));
                    siOrg(i,:) = (siOrg(i,:)-m) .* amp + m;
%                    figure; hold on; plot(signals{k}(i,:)'); plot(siOrg(i,:)'); hold off; title(['node' num2str(i) ' amp=' num2str(amp)]);
%}
                    [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(siOrg);
%%{
                    figure; hold on; plot(si','Color',[0.8 0.8 0.8]); plot(si(i,:)'); hold off; 
%%}
                    for n=1:nMax % traial
%                    parfor n=1:nMax % traial
                        % train DLCM network with amplitude expanded signal
                        netDLCM = initDlcmNetwork(si, f.exSignal, [], f.exControl); 

                        disp(['training 1st ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);
                        netDLCM = trainDlcmNetwork(si, f.exSignal, [], f.exControl, netDLCM, options);

                        % simulate signal from first frame
                        [si2, time] = simulateDlcmNetwork(si, f.exSignal, [], f.exControl, netDLCM);
                        
                        % train DLCM network with simulated signal of amplitude expanded DLCM network
                        netDLCM2 = initDlcmNetwork(si2, f.exSignal, [], f.exControl); 

                        disp(['training 2nd ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);
                        netDLCM2 = trainDlcmNetwork(si2, f.exSignal, [], f.exControl, netDLCM2, options);
                        
                        % calculate DLCM-EC
                        [~, subEC2s{i,a,n}] = calcDlcmEC(netDLCM2, [], f.exControl);
                    end
                    for n=1:nMax
                        Zi2(i, a, n) = subEC2s{i,a,n}(i,1);
                        Zij2(i, a, n,:) = subEC2s{i,a,n}(i,2:end);
                        X(i, a, n) = amp;
                    end
                    save(tmpfName, 'subEC2s', 'Zi2', 'Zij2', 'X');
                end
            end
            save(outfName, 'subEC2s', 'Zi2', 'Zij2', 'X');
            
            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end

        % calc & plot correlation of original Zij vs simulating Zij
        plotCorrelationZiZij([], subEC, [], smSubEC, R, ['sbj' num2str(k)], 'original', 'shifted sim');
        for i=1:R
            for a=1:ampsLen
                for n=1:nMax % traial
                    plotCorrelationZiZij([], subEC, [], subEC2s{i,a,n}, R, ['sbj' num2str(k)], 'original', 'shifted sim');
                end
            end
        end
    end
%}

    % checking signal amplitude change effect for Zi, Zij and ECij'
    global dlcmInitWeights;
    amps = [1, 1.5, 2, 3, 4];
    ampsLen = length(amps);
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);
        smEC = smDLWs(:,:,k);
        smSubEC = smSubDLWs(:,:,k);
        
        outfName = ['results/adsim2-checkRelation5-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            Zi2 = zeros(R, ampsLen, nMax);
            X = zeros(R, ampsLen, nMax);
            Zij2 = zeros(R, ampsLen, nMax, nodeNum);
            subEC2s = cell(R,ampsLen,nMax);
            
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 11;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility

            % training options for DLCM network
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

            for a=4:4 %1:ampsLen
                n = 1;
                amp = amps(a);

                % train DLCM network with amplitude expanded signal
                netDLCM = initDlcmNetwork(signals{k}, f.exSignal, [], f.exControl); 

                nodeLayers = netDLCM.nodeLayers;
                nodeNetwork = cell(nodeNum,1);
                trainInfo = cell(nodeNum,1);
                initWeights = cell(nodeNum,1);

%                for i=1:nodeNum
                parfor i=1:nodeNum
                    siOrg = signals{k};
                    m = nanmean(siOrg(i,:));
                    siOrg(i,:) = (siOrg(i,:)-m) .* amp + m;
                    [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(siOrg);
%{
                    figure; hold on; plot(si','Color',[0.8 0.8 0.8]); plot(si(i,:)'); hold off; 
%}
                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); f.exSignal(:,1:end-1)];
                    filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    disp(['training 1st ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);
                    [nodeNetwork{i}, trainInfo{i}] = trainNetwork(nodeInput, nodeTeach, nodeLayers{i}, options);
                    initWeights{i} = dlcmInitWeights;
                end
                netDLCM.nodeNetwork = nodeNetwork;
                netDLCM.trainInfo = trainInfo;
                netDLCM.initWeights = initWeights;
                netDLCM.trainOptions = options;

                [si2, time] = simulateDlcmNetwork(signals{k}, f.exSignal, [], f.exControl, netDLCM);
                        
                i = 1;
                % train DLCM network with simulated signal of amplitude expanded DLCM network
                disp(['training 2nd ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);

                netDLCM2 = initDlcmNetwork(si2, f.exSignal, [], f.exControl); 
                nodeTeach = si2(i,2:end);
                nodeInput = [si2(:,1:end-1); f.exSignal(:,1:end-1)];
                filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
                
                [trainedNet, ~] = trainNetwork(nodeInput, nodeTeach, netDLCM2.nodeLayers{i}, options);
                        
                % predict DLCM network
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                subEC2s{i,a,n} = predict(trainedNet, Si1);
            end
            for n=1:nMax
                Zi2(i, a, n) = subEC2s{i,a,n}(i,1);
                Zij2(i, a, n,:) = subEC2s{i,a,n}(i,2:end);
                X(i, a, n) = amp;
            end

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
            save(outfName, 'subEC2s', 'Zi2', 'Zij2', 'X');
        end

        % plot correlation of original Zij vs simulating Zij
        plotCorrelationZiZij([], subEC, [], smSubEC, R, ['sbj' num2str(k)], 'original', 'shifted sim');
        for i=1:1 %R
            for a=4:4 %1:ampsLen
                for n=1:1 %nMax % traial
                    plotCorrelationZiZij([], subEC, [], subEC2s{i,a,n}, R, ['sbj' num2str(k)], 'original', 'shifted sim');
                end
            end
        end
    end
end

function checkRelationSubDLWandSignals(signals, DLWs, subDLWs, group, isRaw)
    nodeNum = size(signals{1},1);
    sigLen = size(signals{1},2);
    sbjNum = length(signals);
    R = nodeNum;
    nMax = 10;
    sbjMax = 4;

    % checking signal parallel shift effect for Zi, Zij and ECij'
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);
        
        outfName = ['results/adsim2-checkRelation-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 10;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end

            tmpfName = ['results/adsim2-checkRelation-' group '-' num2str(k) '_tmp.mat'];
            if exist(tmpfName, 'file')
                load(tmpfName);
                if size(X,1) < R
                    X(R,nMax,17) = 0;
                end
                idx = find(X(:, 1, 1)==0);
                rstart =idx(1);
            else
                Zi2 = zeros(R, nMax, 17);
                X = zeros(R, nMax, 17);
                Zij2 = zeros(R, nMax, 17, nodeNum);
                rstart = 1;
            end
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            if isRaw
                siOrg = signals{k};
            else
                [siOrg, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{k});
            end
            % training options for DLCM network
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
                
            for i=rstart:R 
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                for a=0:16
                    dx = (-0.32 + 0.04*a);
                    si = siOrg;
                    si(i,:) = si(i,:) + dx;

                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); f.exSignal(:,1:end-1)];
                    filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    subEC2 = cell(nMax,1);
%                    for n=1:nMax % traial
                    parfor n=1:nMax % traial
                        netDLCM = initDlcmNetwork(si, f.exSignal, [], f.exControl); 

                        disp(['training ' num2str(k) '-' num2str(i) ' dx=' num2str(dx) ' n:' num2str(n)]);
                        [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                        % predict DLCM network
                        subEC2{n} = predict(nodeNetwork, Si1);
                    end
                    for n=1:nMax
                        Zi2(i, n, a+1) = subEC2{n}(1);
                        Zij2(i, n, a+1,:) = subEC2{n}(2:end);
                        X(i, n, a+1) = dx;
                    end
                end
                save(tmpfName, 'Zi2', 'Zij2', 'X');
            end
            save(outfName, 'Zi2', 'Zij2', 'X');

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
        
        % plot result -- Zi2 vs dx
        figure; hold on;
        for i=1:R
            x=X(i,:,:);
            y=Zi2(i,:,:);
            scatter(x(:),y(:),3);
        end
        a=1.5; b=nanmean(Zi2(:,:,9),'all');
        plot([-0.4 0.4], [-0.4*a+b 0.4*a+b],':','Color',[0.2 0.2 0.2]);
        hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zi vs dx']);

        % plot result -- Zij2(1:16) vs dx
        figure; hold on; 
        for i=1:16
            for j=1:R
                x=X(i,:,:);
                y=Zij2(i,:,:,j);
                scatter(x(:),y(:),3);
            end
        end
        plot([-0.4 0.4], [-0.4*a+b 0.4*a+b],':','Color',[0.2 0.2 0.2]);
        hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zij vs dx']);

        % plot result -- Zi - Zij2(1:16) vs dx
        figure; hold on;
        for i=1:16
            for j=1:R
                x=X(i,:,:);
                y=Zi2(i,:,:) - Zij2(i,:,:,j);
                scatter(x(:),y(:),3);
            end
        end
        hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' (Zi - Zij) vs dx']);
    end
    % checking signal amplitude change effect for Zi, Zij and ECij'
%{
    amps = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5, 2, 3, 5, 8];
    ampsLen = length(amps);
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);

        outfName = ['results/adsim2-checkRelation2-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 10;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end
    
            Zi2 = zeros(R, nMax, ampsLen);
            X = zeros(R, nMax, ampsLen);
            Zij2 = zeros(R, nMax, ampsLen, nodeNum);
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            [siOrg, sig, c, maxsi, minsi] = convert2SigmoidSignal(rawSignals{k});

            % training options for DLCM network
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
                
            for i=1:R 
                Zi(i) = subEC(i,1); % original Zi value
                Zij(i,:) = subEC(i,2:end); % original Zij value
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                for a=1:ampsLen
                    si = siOrg;
                    amp = amps(a);
                    m = nanmean(si(i,:));
                    si(i,:) = (si(i,:)-m) .* amp + m;

                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); f.exSignal(:,1:end-1)];
                    filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    subEC2 = cell(nMax,1);
%                    for n=1:nMax % traial
                    parfor n=1:nMax % traial
                        netDLCM = initDlcmNetwork(si, f.exSignal, [], f.exControl); 

                        disp(['training ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);
                        [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                        % predict DLCM network
                        subEC2{n} = predict(nodeNetwork, Si1);
                    end
                    for n=1:nMax
                        Zi2(i, n, a) = subEC2{n}(1);
                        Zij2(i, n, a,:) = subEC2{n}(2:end);
                        X(i, n, a) = amp;
                    end
                end
            end
            save(outfName, 'Zi', 'Zi2', 'Zij2', 'X');

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
        
        % plot result -- Zi2 vs dx
        figure; hold on; 
        for i=1:R
            x=X(i,:,:);
            y=Zi2(i,:,:);
            scatter(x(:),y(:),3);
        end
        hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' Zi vs dx']);

        % plot result -- Zij2(1:64) vs dx
        for i=1:1
            figure; hold on;
            for j=1:64
                x=X(i,:,:);
                y=Zij2(i,:,:,j);
                scatter(x(:),y(:),3);
            end
            hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' node' num2str(i) ' Zij vs dx']);
        end

        % plot result -- Zi - Zij2(1:64) vs dx
        for i=1:1 %R
            figure; hold on;
            for j=1:64
                x=X(i,:,:);
                y=Zi2(i,:,:) - Zij2(i,:,:,j);
                scatter(x(:),y(:),3); 
            end
            Y = repmat(squeeze(Zi2(i,:,:)),[1 1 nodeNum]) - squeeze(Zij2(i,:,:,:));
            ecd1 = Y(:,7,:);
            ecd8 = Y(:,ampsLen,:);
            sd1 = std(ecd1(:),1);
            sd8 = std(ecd8(:),1);
            epr = 8;
            a=0.1*(sd8/sd1)/epr; b=0;
            a2(i) = a;
            plot([0.5 8], [0.5*a+b 8*a+b],':','Color',[0.2 0.2 0.2]);
            plot([0.5 8], [-0.5*a+b -8*a+b],':','Color',[0.2 0.2 0.2]);
            hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' node' num2str(i) ' (Zi - Zij) vs dx']);
        end
    end
%}
    % checking signal amplitude change effect for Zi, Zij and ECij'
    amps = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.2, 1.5, 2, 3, 5, 8];
    ampsLen = length(amps);
    for k=1:sbjMax
        EC = DLWs(:,:,k);
        subEC = subDLWs(:,:,k);

        outfName = ['results/adsim2-checkRelation3-' group '-' num2str(k) '.mat'];
        if exist(outfName, 'file')
            load(outfName);
        else
            % if you want to use parallel processing, set NumProcessors more than 2
            % and change for loop to parfor loop
            NumProcessors = 10;

            if NumProcessors > 1
                try
                    disp('Destroing any existance matlab pool session');
                    parpool('close');
                catch
                    disp('No matlab pool session found');
                end
                parpool(NumProcessors);
            end

            tmpfName = ['results/adsim2-checkRelation3-' group '-' num2str(k) '_tmp.mat'];
            if exist(tmpfName, 'file')
                load(tmpfName);
                if size(X,1) < R
                    X(R,nMax,ampsLen) = 0;
                end
                idx = find(X(:, 1, 3)==0);
                rstart =idx(1);
            else
                Zi2 = zeros(R, nMax, ampsLen);
                X = zeros(R, nMax, ampsLen);
                Zij2 = zeros(R, nMax, ampsLen, nodeNum);
                rstart = 1;
            end
            
            dlcmName = ['results/ad-dlcm-' group '-roi' num2str(nodeNum) '-net' num2str(k) '.mat'];
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            if isRaw
                siOrg = signals{k};
            else
                [siOrg, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{k});
            end

            % training options for DLCM network
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
                
            for i=rstart:R 
                Si1 = ones(nodeNum*2, nodeNum+1);
                Si1(1:nodeNum, 2:end) = ones(nodeNum,nodeNum) - eye(nodeNum);
                filter = repmat(f.exControl(i,:).', 1, size(Si1,2));
                Si1(nodeNum+1:end,:) = Si1(nodeNum+1:end,:) .* filter;

                for a=1:ampsLen
                    si = siOrg;
                    mvsi = movmean(si,5,2);
                    amp = amps(a);
                    si(i,:) = (si(i,:)-mvsi(i,:)) .* amp + mvsi(i,:);
%{
                    figure; hold on; plot(siOrg(i,:)'); plot(mvsi(i,:)'); plot(si(i,:)'); hold off; title(['node' num2str(i) ' amp=' num2str(amp)]);
%}
                    nodeTeach = si(i,2:end);
                    nodeInput = [si(:,1:end-1); f.exSignal(:,1:end-1)];
                    filter = repmat(f.exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;

                    subEC2 = cell(nMax,1);
%                    for n=1:nMax % traial
                    parfor n=1:nMax % traial
                        netDLCM = initDlcmNetwork(si, f.exSignal, [], f.exControl); 

                        disp(['training ' num2str(k) '-' num2str(i) ' amp=' num2str(amp) ' n:' num2str(n)]);
                        [nodeNetwork, trainInfo] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{i}, options);

                        % predict DLCM network
                        subEC2{n} = predict(nodeNetwork, Si1);
                    end
                    for n=1:nMax
                        Zi2(i, n, a) = subEC2{n}(1);
                        Zij2(i, n, a,:) = subEC2{n}(2:end);
                        X(i, n, a) = amp;
                    end
                end
                save(tmpfName, 'Zi2', 'Zij2', 'X');
            end
            save(outfName, 'Zi2', 'Zij2', 'X');

            % shutdown parallel processing
            if NumProcessors > 1
                delete(gcp('nocreate'))
            end
        end
        
        % plot result -- Zi2 vs dx
        figure; hold on; 
        for i=1:R
            x=X(i,:,:);
            y=Zi2(i,:,:);
            scatter(x(:),y(:),3);
        end
        hold off; title(['sbj' num2str(k) ' Zi vs dx']);

        % plot result -- Zij2(1:64) vs dx
        for i=1:4
            figure; hold on;
            for j=1:R
                x=X(i,:,:);
                y=Zij2(i,:,:,j);
                scatter(x(:),y(:),3);
            end
            hold off; title(['sbj' num2str(k) ' node' num2str(i) ' Zij vs dx']);
        end

        % plot result -- Zi - Zij2(1:64) vs dx
%{
        for i=1:R
            figure; hold on;
            for j=1:64
                x=X(i,:,:);
                y=Zi2(i,:,:) - Zij2(i,:,:,j);
                scatter(x(:),y(:),3); 
            end
            hold off; daspect([1 1 1]); title(['sbj' num2str(k) ' node' num2str(i) ' (Zi - Zij) vs dx']);
        end
%}
        x = zeros(R, ampsLen);
        y = zeros(R, ampsLen);
        figure; hold on;
        for i=1:R
            ECds = repmat(squeeze(Zi2(i,:,:)),[1 1 nodeNum]) - squeeze(Zij2(i,:,:,:));
            for a=1:ampsLen
                x(i,a)=X(i,1,a);
                ecd=ECds(:,a,:);
                y(i,a)=std(ecd(:),1);
            end
            scatter(x(i,:),y(i,:), 3, [0.1*mod(i,10) 0.1*ceil(mod(i,100)/10) 0.3+0.2*ceil(i/100)]);
            plot(x(i,:),y(i,:), ':', 'Color', [0.1*mod(i,10) 0.1*ceil(mod(i,100)/10) 0.3+0.2*ceil(i/100)]);
        end
        mx=mean(x,1); my=mean(y,1);
        scatter(mx,my, 7, [0 0 0], 'd');
        plot(mx,my, '--', 'Color', [0 0 0], 'LineWidth',0.5);
        a=0.032; b=0.005;
        plot([1 8], [1*a+b 8*a+b],'-','Color',[1 0.2 0.2]);
        hold off; title(['sbj' num2str(k) ' node' num2str(i) ' std(Zi - Zij) vs dx']);
    end
end

function [smDLWs, bSignals] = retrainDLCMAndECmultiPattern(rawSignals, DLWs, subDLWs, S2, IS2, roiNames, group)
    nodeNum = size(rawSignals{1},1);
    sigLen = size(rawSignals{1},2);
    sbjNum = length(rawSignals);
    nanx = eye(nodeNum);
    nanx(nanx==1) = NaN;

    R = nodeNum;
    JMAX = 7;
    k1 = floor(101/20)+1;
    r1 = zeros(JMAX+1,k1,R);
    r2 = zeros(JMAX+1,k1,R,sbjNum);
    r3 = zeros(JMAX+1,k1,sbjNum);
    h1 = zeros(JMAX+1,k1,nodeNum,nodeNum);
    p1 = zeros(JMAX+1,k1,nodeNum,nodeNum);
    h1c = zeros(JMAX+1,k1);
    r1b = zeros(JMAX+1,k1,R);
    r2b = zeros(JMAX+1,k1,R,sbjNum);
    r3b = zeros(JMAX+1,k1,sbjNum);
    h1b = zeros(JMAX+1,k1,nodeNum,nodeNum);
    p1b = zeros(JMAX+1,k1,nodeNum,nodeNum);
    h1bc = zeros(JMAX+1,k1);

    Zi = repmat(subDLWs(:,1,:),[1 nodeNum 1]);
    Zij = subDLWs(:,2:nodeNum+1,:);
    DLWsR = Zi - Zij;
%    DLWs = abs(DLWsR);

    for ii=0:0
        for exRate=2:2 %0:6
            for k=1:1 %1:20:101
                ECij = Zij - DLWsR * exRate * 0.2;
                teaches = [];
                S3 = [];
                IS3 = [];
                signals = {};
                for i=1:sbjNum
                    [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(rawSignals{i});
                    signals{end+1} = si;
                    teach = [repmat(subDLWs(:,1,i),[1 k 1]) ECij(:,:,i)]; % indivisual part of teaching data
                    teaches(:,:,i) = [teach si(:,2:end)];
                    S3(:,:,i) = [repmat(S2(:,1),[1 k]) S2(:,2:nodeNum+1) si(:,1:end-1)];
                    IS3(:,:,i) = [repmat(IS2(:,1),[1 k]) IS2(:,2:nodeNum+1) rand(nodeNum,sigLen-1)];
                end
                
                name = [group '-' num2str(ii) '-' num2str(exRate) '-' num2str(k) 'ns'];
                bname = [group 'b-' num2str(ii) '-' num2str(exRate) '-' num2str(k) 'ns'];
                ecname = [group '-' num2str(ii) '-' num2str(exRate)  '-' num2str(k) 'ec'];

                % train DLCM network by signals & expanded EC
                [bDLWs, meanbDLWns, stdbDLWns, bSubDLWs] = retrainDLCMAndEC(teaches, S3, IS3, roiNames, name);
                
                k1 = floor(k/20)+1;
                for b=1:R
                    r1b(exRate+1,k1,b) = corr2(squeeze(subDLWs(b,1,:)), squeeze(bSubDLWs(b,1,:)));
%                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' name ' row=' num2str(b)]);
                    for a=1:sbjNum
%                        plotTwoSignalsCorrelation(subDLWs(b,1,a), bSubDLWs(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
%                        plotTwoSignalsCorrelation(subDLWs(b,2:77,a), bSubDLWs(b,2:77,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]);
                        r2b(exRate+1,k1,b,a) = corr2(subDLWs(b,2:1+nodeNum,a), bSubDLWs(b,2:1+nodeNum,a));
                    end; hold off;
                end
%                figure; hold on; plot([0 0.5], [0 0.5],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' ecname ' row=' num2str(b)]);
                for a=1:sbjNum
                    X = DLWs(1:R,1:R,a)+nanx(1:R,1:R);
                    Y = bDLWs(1:R,1:R,a);
%                    plotTwoSignalsCorrelation(X, Y, [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                    r3b(exRate+1,k1,a) = corr2(X(~isnan(X(:))), Y(~isnan(Y(:))));
                end; hold off;
%                calculateAlzWilcoxonTest(subDLWs, bSubDLWs, roiNames, 'ns', name, 'dlw', 1, 'ranksum');
                [h1b(exRate+1,k1,:,:), p1b(exRate+1,k1,:,:), ~] = calculateAlzWilcoxonTest(DLWs, bDLWs, roiNames, 'ec', ecname, 'dlw', 1, 'ranksum', 0);
                h1bc(exRate+1,k1) = length(find(h1b(exRate+1,k1,1:R,:)>0));

                % simulate from first frame with the DLCM network trained by signals & expanded EC
                [bDLWs, bSignals, bSubDLWs] = simulateNodeSignals(signals, roiNames, name, 'dlw', name, 1, [k+nodeNum+1:k+nodeNum+sigLen-1]);
                [smDLs, meanSmDL, ~] = calculateConnectivity(bSignals, roiNames, bname, 'dlcm', 1);
                [smDLWs, meanSmDLW, ~, smSubDLWs] = calculateConnectivity(bSignals, roiNames, bname, 'dlw', 1);

                for b=1:R
                    r1(exRate+1,k1,b) = corr2(squeeze(subDLWs(b,1,:)), squeeze(smSubDLWs(b,1,:)));
%                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' name ' row=' num2str(b)]);
                    for a=1:sbjNum
%                        plotTwoSignalsCorrelation(subDLWs(b,1,a), smSubDLWs(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
%                        plotTwoSignalsCorrelation(subDLWs(b,2:77,a), smSubDLWs(b,2:77,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]);
                        r2(exRate+1,k1,b,a) = corr2(subDLWs(b,2:1+nodeNum,a), smSubDLWs(b,2:1+nodeNum,a));
                    end; hold off;
                end
%                figure; hold on; plot([0 0.5], [0 0.5],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' ecname ' row=' num2str(b)]);
                for a=1:sbjNum
                    X = DLWs(1:R,1:R,a)+nanx(1:R,1:R);
                    Y = smDLWs(1:R,1:R,a);
%                    plotTwoSignalsCorrelation(X, Y, [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                    r3(exRate+1,k1,a) = corr2(X(~isnan(X(:))), Y(~isnan(Y(:))));
                end; hold off;
%                calculateAlzWilcoxonTest(subDLWs, smSubDLWs, roiNames, 'ns', name, 'dlw', 1, 'ranksum');
                [h1(exRate+1,k1,:,:), p1(exRate+1,k1,:,:), ~] = calculateAlzWilcoxonTest(DLWs, smDLWs, roiNames, 'ec', ecname, 'dlw', 1, 'ranksum', 0);
                h1c(exRate+1,k1) = length(find(h1(exRate+1,k1,1:R,:)>0));
            end
        end
    end
    r1m = nanmean(r1,3);
    r2m = nanmean(nanmean(r2,4),3);
    r3m = nanmean(r3,3);
    p1m = nanmean(nanmean(p1(:,:,1:R,:),4),3);
    r1bm = nanmean(r1b,3);
    r2bm = nanmean(nanmean(r2b,4),3);
    r3bm = nanmean(r3b,3);
    p1bm = nanmean(nanmean(p1b(:,:,1:R,:),4),3);
end

function [weights, meanWeights, stdWeights, subweights] = retrainDLCMAndEC(teachSignals, nodeSignals, exSignals, roiNames, group)
    ROWNUM = size(teachSignals,1);
    COLNUM = size(teachSignals,2);
    sbjNum = size(teachSignals,3);
    weights = zeros(ROWNUM, ROWNUM, sbjNum);

    outfName = ['results/adsim2-retrain-' group '-roi' num2str(ROWNUM) '.mat'];
    if exist(outfName, 'file')
        f=load(outfName);
        weights = f.weights;
        subweights = f.subweights;
        meanWeights = nanmean(f.weights, 3);
        stdWeights = nanstd(f.weights, 1, 3);
        return;
    end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 14;

    if NumProcessors > 1
        try
            disp('Destroing any existance matlab pool session');
            parpool('close');
        catch
            disp('No matlab pool session found');
        end
        parpool(NumProcessors);
    end

    % init params
    sigLen = size(nodeSignals,2);
    exControl = eye(ROWNUM);

%    for i=1:sbjNum
    parfor i=1:sbjNum
        dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROWNUM) '-net' num2str(i) '.mat'];
        if exist(dlcmName, 'file')
            f=load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            netDLCM = f.netDLCM;
        else
            if size(nodeSignals,3) > 1
                si = nodeSignals(:,:,i);
                exSignal = exSignals(:,:,i);
            else
                si = nodeSignals;
                exSignal = exSignals;
            end
            % init DLCM network
            netDLCM = initDlcmNetwork(si, exSignal, [], exControl);

            % training DLCM network
            maxEpochs = 1000;
            miniBatchSize = ceil(sigLen / 2);
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
            for j=1:ROWNUM
                nodeTeach = teachSignals(j,1:end,i);
                nodeInput = [si; exSignal];
                if ~isempty(exControl)
                    filter = repmat(exControl(i,:).', 1, size(nodeInput,2));
                    nodeInput(ROWNUM+1:end,:) = nodeInput(ROWNUM+1:end,:) .* filter;
                end
                idx = find(isnan(nodeTeach));
                nodeTeach(:,idx) = [];
                nodeInput(:,idx) = [];
                [netDLCM.nodeNetwork{j}, netDLCM.trainInfo{j}] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{j}, options);
                disp(['virtual alzheimer (' group ') training node ' num2str(i) '-' num2str(j) ' rmse=' num2str(netDLCM.trainInfo{j}.TrainingRMSE(maxEpochs))]);
            end

            parsavedlsm(dlcmName, netDLCM, si, exSignal, exControl, options);
        end

        % recalculate EC
        [weights(:,:,i), subweights(:,:,i)] = calcDlcmEC(netDLCM, [], exControl);
    end
    save(outfName, 'weights', 'roiNames', 'subweights');
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
end

function parsavedlsm(dlcmName, netDLCM, si, exSignal, exControl, options)
    save(dlcmName, 'netDLCM', 'si', 'exSignal', 'exControl', 'options');
end

function out = expandAmplitude(signals, rate)
    m = nanmean(signals(:));
    d = signals - m;
    out = d * rate + m;
end

function out = expandAmplitude2(signals, rate)
    m = nanmean(signals,2);
    d = signals - m;
    out = d * rate + m;
end

function [ECs, simSignals, subECs] = simulateNodeSignals(signals, roiNames, group, algorithm, orgGroup, isRaw, inSiRange)
    if nargin < 7
        inSiRange = 0;
    end
    if nargin < 6
        isRaw = 0;
    end
    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    outfName = ['results/adsim2-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end

    % if you want to use parallel processing, set NumProcessors more than 2
    % and change for loop to parfor loop
    NumProcessors = 11;

    if NumProcessors > 1
        try
            disp('Destroing any existance matlab pool session');
            parpool('close');
        catch
            disp('No matlab pool session found');
        end
        parpool(NumProcessors);
    end
        
    ECs = zeros(ROINUM, ROINUM, sbjNum);
    subECs = zeros(ROINUM, ROINUM+1, sbjNum);
    simSignals = cell(1, sbjNum);
    parfor i=1:sbjNum    % for parallel processing
%    for i=1:sbjNum
        switch(algorithm)
        case {'dlw','dlwrc'}
            if strcmp(algorithm, 'dlw')
                dlcmName = ['results/ad-dlcm-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            else
                dlcmName = ['results/ad-dlcmrc-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
            end
            f = load(dlcmName);
            if isfield(f,'inSignal'), f.exSignal = f.inSignal; end % for compatibility
            if isfield(f,'inControl'), f.exControl = f.inControl; end % for compatibility
            if isRaw
                si = signals{i};
                sig=0; c=0; maxsi=0; minsi=0;
            else
                [si, sig, c, maxsi, minsi] = convert2SigmoidSignal(signals{i});
            end
            if inSiRange > 0
                exSignal = f.exSignal(:,inSiRange);
            else
                exSignal = f.exSignal;
            end
            [Y, time] = simulateDlcmNetwork(si, exSignal, [], f.exControl, f.netDLCM);
            [ec, subECs(:,:,i)] = calcDlcmEC(f.netDLCM, [], f.exControl);
        end
        ECs(:,:,i) = ec;
        simSignals{i} = Y;
    end
    save(outfName, 'ECs', 'simSignals', 'roiNames', 'subECs');

    % shutdown parallel processing
    if NumProcessors > 1
        delete(gcp('nocreate'))
    end
end


function [utestH, utestP, utestP2] = calculateAlzWilcoxonTest(control, target, roiNames, controlGroup, targetGroup, algorithm, force, testname, showfig)
    if nargin < 9
        showfig = 1;
    end
    if nargin < 8
        testname = 'ranksum';
    end
    if nargin < 7
        force = 0;
    end
    % constant value
    ROWNUM = size(control,1);
    COLNUM = size(control,2);

    outfName = ['results/adsim-' algorithm '-' controlGroup '_' targetGroup '-roi' num2str(ROWNUM) '-utest.mat'];
    if exist(outfName, 'file') && force == 0
        load(outfName);
    else
        utestH = nan(ROWNUM, COLNUM);
        utestP = nan(ROWNUM, COLNUM);
        utestP2 = nan(ROWNUM, COLNUM);
        for i=1:ROWNUM
            for j=1:COLNUM
                if isnan(control(i,j,1)), continue; end
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
                switch(testname)
                    case 'ranksum'
                        [p, h] = ranksum(x2,y2);
                    case 'ttest2'
                        [h, p] = ttest2(x2,y2);
                end
                %[p, h] = signrank(x2,y2);
                utestH(i,j) = h;
                utestP(i,j) = p;
                if h > 0 && nanmean(x) > nanmean(y)
                    utestP2(i,j) = p;
                end
            end
        end
        if force == 0
            save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames');
        end
    end
    % counting by source region and target region
    countSource = nansum(utestH,1);
    countTarget = nansum(utestH,2);
    save(outfName, 'utestH', 'utestP', 'utestP2', 'roiNames', 'countSource', 'countTarget');
    if showfig == 0, return; end
    
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
