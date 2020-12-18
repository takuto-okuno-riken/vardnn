% this function works after analyzeAlzheimerDLCM.m

function simulateAlzheimerDLCM
    % CONN fmri data base path :
    base = '../fmri/';

    % CONN output path
    pathesCN = {'ADNI2_65-78_F_CN_nii', 'ADNI2_65-78_M_CN_nii'};
    pathesAD = {'ADNI2_65-75_F_AD_nii', 'ADNI2_65-75_M_AD_nii'};

    % load each type signals
    [cnSignals, roiNames] = connData2signalsFile(base, pathesCN, 'cn');
    [adSignals] = connData2signalsFile(base, pathesAD, 'ad');

    % calculate node-signals from (1,...,1), (0,1...1)->(1...1,0),
    % (0,...,0), (1,0...0)->(0...0,1), (0.5,0...0)->(0...0,0.5)
%    [cnDLWs, cnDLWnss, meanCnDLWns, stdCnDLWns, cnInSignals, cnInControls] = calculateDistributions(cnSignals, roiNames, 'cn', 'dlw');
%    [adDLWs, adDLWnss, meanAdDLWns, stdAdDLWns, ~, ~] = calculateDistributions(adSignals, roiNames, 'ad', 'dlw');
    [cnDLWs, cnDLWnss, meanCnDLWns, stdCnDLWns, cnInSignals, cnInControls, cnS2, cnIS2] = calculateDistributions2(cnSignals, roiNames, 'cn', 'dlw', 'cn');
    [adDLWs, adDLWnss, meanAdDLWns, stdAdDLWns, ~, ~, ~, ~] = calculateDistributions2(adSignals, roiNames, 'ad', 'dlw', 'ad');

    % transform healthy node signals to ad's distribution (type 1)
    ROINUM = size(cnDLWs, 1);
    cnSbjNum = size(cnDLWs,3);
    adSbjNum = size(adDLWs,3);
    cnZi = repmat(cnDLWnss(:,1,:),[1 ROINUM 1]);
    cnZij = cnDLWnss(:,2:ROINUM+1,:);
    cnDLWsR = (cnZi - cnZij);   % non-abs EC of AD

    adZi = repmat(adDLWnss(:,1,:),[1 ROINUM 1]);
    adZij = adDLWnss(:,2:ROINUM+1,:);
    adDLWsR = (adZi - adZij);   % non-abs EC of AD

%    cnDLWsRstd = nanstd(cnDLWsR,1,3);
%    adDLWsRstd = nanstd(adDLWsR,1,3);
%    sigWeight = adDLWsRstd ./ cnDLWsRstd;

%    adPMask = (adDLWsR >= 0);
%    adMMask = (adDLWsR < 0);
%    meanAdZi = mean(adZi,3);
%    meanAdZij = nanmean(adZij,3);
    
    meanCnDLWns3 = repmat(meanCnDLWns,[1 1 cnSbjNum]);
    stdCnDLWns3 = repmat(stdCnDLWns,[1 1 cnSbjNum]);
    sigCnDLWns = (cnDLWnss - meanCnDLWns3) ./ stdCnDLWns3;
    meanAdDLWns3 = repmat(meanAdDLWns,[1 1 cnSbjNum]);
    stdAdDLWns3 = repmat(stdAdDLWns,[1 1 cnSbjNum]);
    vadDLWstd = sigCnDLWns .* stdAdDLWns3;
    vadDLWnss = meanAdDLWns3 + vadDLWstd;

    nanx = eye(ROINUM);
    nanx(nanx==1) = NaN;

    % --------------------------------------------------------------------------------------------------------------
    % calculate virtual AD ECcnDLWnss (type 1 : EC, teach-signals)
%    vadDLWs = adPMask .* repmat((meanAdZi-meanAdZij),[1 1 adSbjNum]) + adMMask .* repmat((meanAdZij-meanAdZi),[1 1 adSbjNum]);
%    
    vadZi = repmat(vadDLWnss(:,1,:),[1 ROINUM 1]);
    vadZij = vadDLWnss(:,2:ROINUM+1,:);

%    this needs to calc vadDLWsRstd, but it is redundant.
    vadDLWsR = (vadZi - vadZij);   % non-abs EC of virtual AD
%    vadDLWsRstd = nanstd(vadDLWsR,1,3);
%    sigWeight = adDLWsRstd ./ vadDLWsRstd; %

    vadCov = nan(ROINUM,ROINUM);
    for i=1:ROINUM
        for j=1:ROINUM
            X = squeeze(vadZi(i,j,:));
            Y = squeeze(vadZij(i,j,:));
            C = cov(X,Y,1);
            vadCov(i,j) = C(1,2);
        end
    end
    b4ac = 4 .* vadCov .* vadCov - 4 .* nanvar(vadZij,1,3) .* (nanvar(vadZi,1,3) - nanvar(adDLWsR,1,3));
    b4ac(b4ac<0) = 0;
    b4ac = sqrt(b4ac);
    sigWeight = (2 .* vadCov + b4ac) ./ (2 .* nanvar(vadZij,1,3));

    vadDLWnss(:,2:ROINUM+1,:) = meanAdDLWns3(:,2:ROINUM+1,:) +  vadDLWstd(:,2:ROINUM+1,:) .* repmat(sigWeight,[1 1 cnSbjNum]);
    % this expression is better, but Zi from input (1...1) can take only one value
%   vadZi = repmat(meanAdDLWns3(:,1,:),[1 ROINUM 1]) + repmat(vadDLWstd(:,1,:),[1 ROINUM 1]) .* repmat(sigWeight,[1 1 cnSbjNum]);
    % one value of Zi. mean of sigWeight doesn't affect so much
    vadDLWnss(:,1,:) = meanAdDLWns3(:,1,:) + vadDLWstd(:,1,:);% .* repmat(mean(sigWeight,2),[1 1 cnSbjNum]);
    vadZi = repmat(meanAdDLWns3(:,1,:),[1 ROINUM 1]) + repmat(vadDLWstd(:,1,:),[1 ROINUM 1]);

    vadZij = vadDLWnss(:,2:ROINUM+1,:);
    vadDLWsR = (vadZi - vadZij);   % non-abs EC of virtual AD
    vadDLWs = abs(vadDLWsR); % EC of virtual AD

%    vadDLWs = vadDLWsR;
%    vadDLWs(vadDLWs<0) = 0;
%{
i=5; j=2;
a = [];
a(:,1) = squeeze(cnDLWs(i,j,:));
a(:,2) = squeeze(cnZi(i,j,:));
a(:,3) = squeeze(cnZij(i,j,:));

b = [];
b(:,1) = squeeze(adDLWs(i,j,:));
b(:,2) = squeeze(adZi(i,j,:));
b(:,3) = squeeze(adZij(i,j,:));
%}
    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 2 : EC, net)
    [vad2DLWs, meanVad2DLWns, stdVad2DLWns] = retrainDLCMAndEC(vadDLWnss, cnS2, cnIS2, roiNames, 'vadns');
    [vad2bDLWs, vad2DLWnss] = calculateNodeSignals(cnSignals, cnS2, cnIS2, roiNames, 'vadns', 'dlw');
    vad2Zi = repmat(vad2DLWnss(:,1,:),[1 ROINUM 1]);
    vad2Zij = vad2DLWnss(:,2:ROINUM+1,:);
    vad2DLWsR = (vad2Zi - vad2Zij);
%%{
    p = 0;
    for b=1:2
        figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: vad-vad2 row=' num2str(b)]);
        for a=1:cnSbjNum
            plotTwoSignalsCorrelation(vadDLWnss(b,1,a), vad2DLWnss(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
            plotTwoSignalsCorrelation(vadDLWnss(b,p+1:p+66,a), vad2DLWnss(b,p+1:p+66,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]); 
        end; hold off;
    end
    figure; hold on; plot([0 0.2], [0 0.2],':','Color',[0.5 0.5 0.5]); title(['ec corr: vad-vad2 row=' num2str(b)]);
    for a=1:cnSbjNum
        plotTwoSignalsCorrelation(vadDLWs(1:4,1:4,a)+nanx(1:4,1:4), vad2bDLWs(1:4,1:4,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
    end; hold off;
%%}
    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 12 : EC, net) (optimise for DLCM training)
%    [r1m, r2m, r3m, h1c, p1m] = retrainDLCMAndECmultiPattern(cnSignals, adDLWs, vadDLWs, vadDLWnss, vadZij, vadDLWsR, cnS2, cnIS2, roiNames, 'vad18');

    % --------------------------------------------------------------------------------------------------------------
    % transform healthy node signals to ad's distribution (type 13 : EC, teach-signals)
    % first generate vad Zi, then calculate Zij from non-abs EC of AD
    meanAdDLWsR = mean(adDLWsR,3);
    stdAdDLWsR = nanstd(adDLWsR,1,3);
    meanCnDLWsR = mean(cnDLWsR,3);
    stdCnDLWsR = nanstd(cnDLWsR,1,3);
    sigCnDLWsR = (cnDLWsR - repmat(meanCnDLWsR, [1 1 cnSbjNum])) ./ repmat(stdCnDLWsR, [1 1 cnSbjNum]);

%    vad19DLWnss = meanAdDLWns3 + sigCnDLWns .* stdAdDLWns3;
%    vad19Zi = repmat(vad19DLWnss(:,1,:),[1 ROINUM 1]);
%    vad19Zij = vad19DLWnss(:,2:ROINUM+1,:);
%    vad19DLWs = abs(vad19Zi - vad19Zij);

    vad19DLWnss = meanAdDLWns3 + sigCnDLWns .* stdAdDLWns3;
    vad19Zi = repmat(vad19DLWnss(:,1,:),[1 ROINUM 1]);
    vad19Zij = vad19Zi - (repmat(meanAdDLWsR, [1 1 cnSbjNum]) + repmat(stdAdDLWsR, [1 1 cnSbjNum]) .* sigCnDLWsR .* 1.0);
    
    vad19DLWnss(:,2:ROINUM+1,:) = vad19Zij;
    vad19DLWsR = vad19Zi - vad19Zij;
    vad19DLWs = abs(vad19DLWsR);

%    [vad19H, vad19P, ~] = calculateAlzWilcoxonTest(adDLWs, vad19DLWs, roiNames, 'adec', 'vad19ec', 'dlw', 1, 'ranksum');
%    [vadH, vadP, ~] = calculateAlzWilcoxonTest(adDLWs, vadDLWs, roiNames, 'adec', 'vadec', 'dlw', 1, 'ranksum');

    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 14 : EC, net) (optimise for DLCM training)
%    [r1m, r2m, r3m, h1c, p1m] = retrainDLCMAndECmultiPattern(cnSignals, adDLWs, vad19DLWs, vad19DLWnss, vad19Zij, vad19DLWsR, cnS2, cnIS2, roiNames, 'vad21');

    % --------------------------------------------------------------------------------------------------------------
    % transform healthy node signals to ad's distribution (type 15 : EC, teach-signals)
    % first generate vad Zi, then calculate Zij from non-abs EC of AD
    % calculate node-signals other pattern (1,...,1),(0.75,...,0.75),..,(0,...,0)
    % , (1,0,1..0),(0,1,0..1),..,(0.25,0,0.25..0),(0,0.25,0..0.25)

    [cn3DLWs, cn3DLWnss, meanCn3DLWns, stdCn3DLWns, cn3InSignals, cn3InControls, cnS3, cnIS3] = calculateDistributions3(cnSignals, roiNames, 'cn3', 'dlw', 'cn');
    [ad3DLWs, ad3DLWnss, meanAd3DLWns, stdAd3DLWns, ~, ~, ~, ~] = calculateDistributions3(adSignals, roiNames, 'ad3', 'dlw', 'ad');
    meanCn3DLWns3 = repmat(meanCn3DLWns,[1 1 cnSbjNum]);
    stdCn3DLWns3 = repmat(stdCn3DLWns,[1 1 cnSbjNum]);
    sigCn3DLWns = (cn3DLWnss - meanCn3DLWns3) ./ stdCn3DLWns3;
    meanAd3DLWns3 = repmat(meanAd3DLWns,[1 1 cnSbjNum]);
    stdAd3DLWns3 = repmat(stdAd3DLWns,[1 1 cnSbjNum]);

    cn3Zi = repmat(cn3DLWnss(:,ROINUM+1,:),[1 ROINUM 1]);
    cn3Zij = cn3DLWnss(:,1:ROINUM,:);
    cn3DLWsR = (cn3Zi - cn3Zij);   % non-abs EC of AD
    ad3Zi = repmat(ad3DLWnss(:,ROINUM+1,:),[1 ROINUM 1]);
    ad3Zij = ad3DLWnss(:,1:ROINUM,:);
    ad3DLWsR = (ad3Zi - ad3Zij);   % non-abs EC of AD

    meanAd3DLWsR = mean(ad3DLWsR,3);
    stdAd3DLWsR = nanstd(ad3DLWsR,1,3);
    meanCn3DLWsR = mean(cn3DLWsR,3);
    stdCn3DLWsR = nanstd(cn3DLWsR,1,3);
    sigCn3DLWsR = (cn3DLWsR - repmat(meanCn3DLWsR, [1 1 cnSbjNum])) ./ repmat(stdCn3DLWsR, [1 1 cnSbjNum]);

    vad22DLWnss = meanAd3DLWns3 + sigCn3DLWns .* stdAd3DLWns3;
    vad22Zi = repmat(vad22DLWnss(:,ROINUM+1,:),[1 ROINUM 1]);
    vad22Zij = vad22Zi - (repmat(meanAd3DLWsR, [1 1 cnSbjNum]) + repmat(stdAd3DLWsR, [1 1 cnSbjNum]) .* sigCn3DLWsR .* 1.0);
    
    vad22DLWnss(:,1:ROINUM,:) = vad22Zij;
    vad22DLWsR = vad22Zi - vad22Zij;
    vad22DLWs = abs(vad22DLWsR);

%    [vad22H, vad22P, ~] = calculateAlzWilcoxonTest(adDLWs, vad22DLWs, roiNames, 'adec', 'vad22ec', 'dlw', 1, 'ranksum');
%    [vadH, vadP, ~] = calculateAlzWilcoxonTest(adDLWs, vadDLWs, roiNames, 'adec', 'vadec', 'dlw', 1, 'ranksum');

    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 16 : EC, net) (optimise for DLCM training)
    [r1m, r2m, r3m, h1c, p1m] = retrainDLCMAndECmultiPattern(cnSignals, adDLWs, vad22DLWs, vad22DLWnss, vad22Zij, vad22DLWsR, cnS3, cnIS3, roiNames, 'vad24');

    % --------------------------------------------------------------------------------------------------------------
    % transform healthy node signals to ad's distribution (type 5 : EC, teach-signals)
    % first generate vad Zi, then calculate Zij from ad EC (random sigma)
    meanAdDLWs = mean(adDLWs,3);
    stdAdDLWs = nanstd(adDLWs,1,3);

    outfName = ['results/adsim-dlw-vad7ns-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        vad7DLWnss = vadDLWnss;
        vad7Zi = repmat(vadDLWnss(:,1,:),[1 ROINUM 1]);
    %    vad7Zij = vad7Zi - repmat(meanAdDLWs, [1 1 cnSbjNum]) + repmat(stdAdDLWs, [1 1 cnSbjNum]) .* randn([ROINUM ROINUM cnSbjNum]) .* 0.7;
        vad7Zij = vad7Zi - repmat(meanAdDLWs, [1 1 cnSbjNum]) + repmat(stdAdDLWs, [1 1 cnSbjNum]) .* (rand([ROINUM ROINUM cnSbjNum])-0.5) .* 2;

        vad7DLWnss(:,2:ROINUM+1,:) = vad7Zij;
        vad7DLWsR = (vad7Zi - vad7Zij);
        vad7DLWs = abs(vad7DLWsR);
        save(outfName, 'vad7DLWnss', 'vad7Zi', 'vad7Zij', 'vad7DLWsR', 'vad7DLWs');
    end
    
    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 6 : EC, net)
    [vad8DLWs, meanVad8DLWns, stdVad8DLWns] = retrainDLCMAndEC(vad7DLWnss, cnS2, cnIS2, roiNames, 'vad8ns');

    % --------------------------------------------------------------------------------------------------------------
    % transform healthy node signals to ad's distribution (type 7 : EC, teach-signals)
    % first generate vad Zi, then calculate Zij from ad EC
    meanCnDLWs = mean(cnDLWs,3);
    stdCnDLWs = nanstd(cnDLWs,1,3);
    sigCnDLW = (cnDLWs - repmat(meanCnDLWs, [1 1 cnSbjNum])) ./ repmat(stdCnDLWs, [1 1 cnSbjNum]);
    
    vad9DLWnss = vadDLWnss;
    vad9Zi = repmat(vadDLWnss(:,1,:),[1 ROINUM 1]);
    vad9Zij = vad9Zi - (repmat(meanAdDLWs, [1 1 cnSbjNum]) + repmat(stdAdDLWs, [1 1 cnSbjNum]) .* sigCnDLW .* 0.8);

    vad9DLWnss(:,2:ROINUM+1,:) = vad9Zij;
    vad9DLWs = abs(vad9Zi - vad9Zij);
%    calculateAlzWilcoxonTest(adDLWs, vad9DLWs, roiNames, 'adec', 'vad9ec', 'dlw', 1, 'ttest2');
%    calculateAlzWilcoxonTest(adDLWs, vad9DLWs, roiNames, 'adec', 'vad9ec', 'dlw', 1, 'ranksum');

    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 8 : EC, net)
    [vad10DLWs, meanVad10DLWns, stdVad10DLWns] = retrainDLCMAndEC(vad9DLWnss, cnS2, cnIS2, roiNames, 'vad10ns');
    [vad10bDLWs, vad10DLWnss] = calculateNodeSignals(cnSignals, cnS2, cnIS2, roiNames, 'vad10ns', 'dlw');
%    figure; plotEC(mean(adDLWs,3),'ad DLW',0);
%    figure; plotEC(mean(vad10DLWs,3),'vad10 DLW',0);
%    [advad10DLWsUt, advad10DLWsUtP, advad10DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad10DLWs, roiNames, 'adec', 'vad10ec', 'dlw', 1, 'ranksum');
%    sigAdDLWs = sigmaEC(mean(adDLWs,3));
%    sigVad10DLWs = sigmaEC(mean(vad10DLWs,3));
%    figure; advadDLWr2 = plotTwoSignalsCorrelation(sigAdDLWs, sigVad10DLWs);
%    sigAdDLWs = sigmaEC(adDLWs);
%    sigCnDLWs = sigmaEC(cnDLWs);
%    sigVad10DLWs = sigmaEC(vad10DLWs);
%    [cnvad10DLWsUt, cnvad10DLWsUtP, cnvad10DLWsUtP2] = calculateAlzWilcoxonTest(sigCnDLWs, sigVad10DLWs, roiNames, 'cnec', 'vad10ec', 'dlw', 1, 'ranksum');
%    [advad10DLWsUt, advad10DLWsUtP, advad10DLWsUtP2] = calculateAlzWilcoxonTest(sigAdDLWs, sigVad10DLWs, roiNames, 'adec', 'vad10ec', 'dlw', 1, 'ranksum');
%    calculateAlzWilcoxonTest(vad9DLWs, vad10DLWs, roiNames, 'vad9ec', 'vad10ec', 'dlw', 1, 'ranksum');
%    calculateAlzWilcoxonTest(vad9DLWnss, vad10DLWnss, roiNames, 'vad9ns', 'vad10ns', 'dlw', 1, 'ranksum');

    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 9 : EC, net) (optimise for DLCM training)
    vad11DLWnss = [repmat(vad9DLWnss(:,1,:),[1 160 1]) vad9DLWnss(:,2:end,:)];
    cnS11 = [repmat(cnS2(:,1,:),[1 160 1]) cnS2(:,2:end,:)];
    cnIS11 = [repmat(cnIS2(:,1,:),[1 160 1]) cnS2(:,2:end,:)];

    [vad12DLWs, meanVad12DLWns, stdVad12DLWns] = retrainDLCMAndEC(vad11DLWnss, cnS11, cnIS11, roiNames, 'vad12ns');
    [vad12bDLWs, vad12DLWnss] = calculateNodeSignals(cnSignals, cnS11, cnIS11, roiNames, 'vad12ns', 'dlw');

%    calculateAlzWilcoxonTest(vad11DLWnss, vad12DLWnss, roiNames, 'vad11ns', 'vad12ns', 'dlw', 1, 'ranksum');
%    calculateAlzWilcoxonTest(vad9DLWs, vad12bDLWs, roiNames, 'vad9ec', 'vad12ec', 'dlw', 1, 'ranksum');
%    calculateAlzWilcoxonTest(vad9DLWs, vad12bDLWs, roiNames, 'vad9ec', 'vad12ec', 'dlw', 1, 'ttest2');

    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 10 : EC, net) (optimise for DLCM training)
    vad13DLWnss = vadDLWnss;
    vad13Zi = repmat(vadDLWnss(:,1,:),[1 ROINUM 1]);
    vad19Zij = vad13Zi - (repmat(meanAdDLWs, [1 1 cnSbjNum]) + repmat(stdAdDLWs, [1 1 cnSbjNum]) .* sigCnDLW .* 0.8);
%{
    for i=0:5
        for j=0:5
            for k=0:3
                p = ceil((ROINUM+1)*(1 + k*0.5));
                vad13DLWnss(:,2:ROINUM+1,:) = vad13Zij - 0.01*i;
                vad13DLWs = abs(vad13Zi - vad13Zij);
                vad13DLWnss(:,1,:) = vadDLWnss(:,1,:) + 0.01*j;
                vad13DLWnss = [repmat(vad13DLWnss(:,1,:),[1 p 1]) vad13DLWnss(:,2:ROINUM+1,:)];
                cnS13 = [repmat(cnS2(:,1,:),[1 p 1]) cnS2(:,2:ROINUM+1,:)];
                cnIS13 = [repmat(cnIS2(:,1,:),[1 p 1]) cnS2(:,2:ROINUM+1,:)];

                vad13name = ['vad13-' num2str(i) '-' num2str(j) '-' num2str(k) 'ns'];
                vad14name = ['vad14-' num2str(i) '-' num2str(j) '-' num2str(k) 'ns'];
                vad14ecname = ['vad14-' num2str(i) '-' num2str(j) '-' num2str(k) 'ec'];
                [vad14DLWs, meanVad14DLWns, stdVad14DLWns] = retrainDLCMAndEC(vad13DLWnss, cnS13, cnIS13, roiNames, vad14name);
                [vad14bDLWs, vad14DLWnss] = calculateNodeSignals(cnSignals, cnS13, cnIS13, roiNames, vad14name, 'dlw');
                
                vad11bDLWnss = [repmat(vad9DLWnss(:,1,:),[1 p 1]) vad9DLWnss(:,2:ROINUM+1,:)];
%                calculateAlzWilcoxonTest(vad11bDLWnss, vad14DLWnss, roiNames, 'vad11ns', vad14name, 'dlw', 1, 'ranksum');
%                calculateAlzWilcoxonTest(vad13DLWnss, vad14DLWnss, roiNames, vad13name, vad14name, 'dlw', 1, 'ranksum');
%                calculateAlzWilcoxonTest(adDLWs, vad14bDLWs, roiNames, 'adec', vad14ecname, 'dlw', 1, 'ranksum');
                for b=1:2
                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' vad14name ' row=' num2str(b)]);
                    for a=1:cnSbjNum
                        plotTwoSignalsCorrelation(vad11bDLWnss(b,1,a), vad14DLWnss(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
                        nss14 = vad14DLWnss(b,p+1:p+40,a); nss14(b) = NaN;
                        plotTwoSignalsCorrelation(vad11bDLWnss(b,p+1:p+40,a), nss14, [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]); 
                    end; hold off;
                end
                figure; hold on; plot([0 0.2], [0 0.2],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' vad14ecname ' row=' num2str(b)]);
                for a=1:cnSbjNum
                    plotTwoSignalsCorrelation(vad9DLWs(1:4,1:4,a), vad14bDLWs(1:4,1:4,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                end; hold off;
            end
        end
    end
%}
    % --------------------------------------------------------------------------------------------------------------
    % re-training DLCM network (type 11 : EC, net) (optimise for DLCM training)
%{
    for i=2:2
        for j=25:25
            for k=2:2
                p = ceil((ROINUM+1)*(k*0.2));
                vad15DLWs = abs(vad13Zi - vad13Zij);
                vad15DLWnss(:,2:ROINUM+1,:) = vad13Zij - vad15DLWs * i * 0.5;
                vad15DLWnss(:,1,:) = vadDLWnss(:,1,:) + nanmean(vad15DLWs,2) * j * 0.2;
                vad15DLWnss = [repmat(vad15DLWnss(:,1,:),[1 p 1]) vad15DLWnss(:,2:ROINUM+1,:)];
                cnS15 = [repmat(cnS2(:,1,:),[1 p 1]) cnS2(:,2:ROINUM+1,:)];
                cnIS15 = [repmat(cnIS2(:,1,:),[1 p 1]) cnS2(:,2:ROINUM+1,:)];

                vad16name = ['vad16-' num2str(i) '-' num2str(j) '-' num2str(k) 'ns'];
                vad16ecname = ['vad16-' num2str(i) '-' num2str(j) '-' num2str(k) 'ec'];
                [vad16DLWs, meanVad16DLWns, stdVad16DLWns] = retrainDLCMAndEC(vad15DLWnss, cnS15, cnIS15, roiNames, vad16name);
                [vad16bDLWs, vad16DLWnss] = calculateNodeSignals(cnSignals, cnS15, cnIS15, roiNames, vad16name, 'dlw');

                vad11bDLWnss = [repmat(vad9DLWnss(:,1,:),[1 p 1]) vad9DLWnss(:,2:ROINUM+1,:)];
    %                calculateAlzWilcoxonTest(vad11bDLWnss, vad14DLWnss, roiNames, 'vad11ns', vad14name, 'dlw', 1, 'ranksum');
    %                calculateAlzWilcoxonTest(vad13DLWnss, vad14DLWnss, roiNames, vad13name, vad14name, 'dlw', 1, 'ranksum');
    %                calculateAlzWilcoxonTest(adDLWs, vad14bDLWs, roiNames, 'adec', vad14ecname, 'dlw', 1, 'ranksum');
                for b=1:2
                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' vad16name ' row=' num2str(b)]);
                    for a=1:cnSbjNum
                        plotTwoSignalsCorrelation(vad11bDLWnss(b,1,a), vad16DLWnss(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
                        nss16 = vad16DLWnss(b,p+1:p+40,a); nss16(b) = NaN;
                        plotTwoSignalsCorrelation(vad11bDLWnss(b,p+1:p+40,a), nss16, [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]); 
                    end; hold off;
                end
                figure; hold on; plot([0 0.2], [0 0.2],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' vad16ecname ' row=' num2str(b)]);
                for a=1:cnSbjNum
                    plotTwoSignalsCorrelation(vad9DLWs(1:4,1:4,a), vad16bDLWs(1:4,1:4,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                end; hold off;
            end
        end
    end
%}
    % --------------------------------------------------------------------------------------------------------------
    % generate virtual ad signals (type 3 : EC, BOLD-signals, net)
    % transform cn signals to vad signals linearly (based on EC rate)
%{
    vadSignals = calculateVirtualADSignals3(cnSignals, roiNames, cnDLWs, adDLWs, 'dlw');
    [vad3DLs, ~, ~] = calculateConnectivity(vadSignals, roiNames, 'vad3', 'dlcm');
    [vad3DLWs, ~, ~] = calculateConnectivity(vadSignals, roiNames, 'vad3', 'dlw');
    [~, vad3DLWnss, meanVad3DLWns, stdVad3DLWns, ~, ~, ~, ~] = calculateDistributions2(vadSignals, roiNames, 'vad3', 'dlw', 'vad3');

    % generate virtual ad signals (type 4 : EC, BOLD-signals, net)
    % input cn-signals to AD-nets and take mean of 32 outputs
    [vad4Signals, vad4DLWs, vad4DLWnss] = calculateVirtualADSignals4(cnSignals, adSignals, roiNames, cnInSignals, cnInControls, 'vad4');
    [vad5Signals, vad5DLWs, vad5DLWnss] = calculateVirtualADSignals4(vad4Signals, adSignals, roiNames, cnInSignals, cnInControls, 'vad5');
%    [vad6Signals, vad6DLWs, vad6DLWnss] = calculateVirtualADSignals4(vad5Signals, adSignals, roiNames, cnInSignals, cnInControls, 'vad6');
%}
    % --------------------------------------------------------------------------------------------------------------
    % calculate node-signals other pattern (1,...,1),(0.75,...,0.75),..,(0,...,0)
    % , (1,0,1..0),(0,1,0..1),..,(0.25,0,0.25..0),(0,0.25,0..0.25)
%{
    [cn3DLWs, cn3DLWnss, meanCn3DLWns, stdCn3DLWns, cn3InSignals, cn3InControls, cnS3, cnIS3] = calculateDistributions3(cnSignals, roiNames, 'cn3', 'dlw', 'cn');
    [ad3DLWs, ad3DLWnss, meanAd3DLWns, stdAd3DLWns, ~, ~, ~, ~] = calculateDistributions3(adSignals, roiNames, 'ad3', 'dlw', 'ad');
    meanCn3DLWns3 = repmat(meanCn3DLWns,[1 1 cnSbjNum]);
    stdCn3DLWns3 = repmat(stdCn3DLWns,[1 1 cnSbjNum]);
    sigCn3DLWns = (cn3DLWnss - meanCn3DLWns3) ./ stdCn3DLWns3;
    meanAd3DLWns3 = repmat(meanAd3DLWns,[1 1 cnSbjNum]);
    stdAd3DLWns3 = repmat(stdAd3DLWns,[1 1 cnSbjNum]);
    vad15DLWstd = sigCn3DLWns .* stdAd3DLWns3;
    vad15DLWnss = meanAd3DLWns3 + vad15DLWstd;

    [vad16DLWs, meanVad16DLWns, stdVad16DLWns] = retrainDLCMAndEC(vad15DLWnss(:,ROINUM+1:end,:), cnS3(:,ROINUM+1:end,:), cnIS3(:,ROINUM+1:end,:), roiNames, 'vad15ns');
    [vad16bDLWs, vad16DLWnss] = calculateNodeSignals(cnSignals, cnS3, cnIS3, roiNames, 'vad15ns', 'dlw');
    
%    calculateAlzWilcoxonTest(vad15DLWnss, vad16DLWnss, roiNames, 'vad15ns', 'vad16ns', 'dlw', 1, 'ranksum');
%    calculateAlzWilcoxonTest(adDLWs, vad16bDLWs, roiNames, 'adec', 'vad16ec', 'dlw', 1, 'ranksum');

%    ad3Zi = repmat(ad3DLWnss(:,ROINUM+1,:),[1 ROINUM 1]);
%    ad3Zij = ad3DLWnss(:,1:ROINUM,:);
%    ad3DLWsR = (ad3Zi - ad3Zij);   % non-abs EC of AD
%    vad16Zij = vad16DLWnss(:,1:ROINUM,:);
%    vad16Zi = repmat(vad16DLWnss(:,ROINUM+1,:),[1 ROINUM 1]);
%    vad16DLWsR = (vad16Zi - vad16Zij);   % non-abs EC of AD
%    calculateAlzWilcoxonTest(ad3DLWsR, vad16DLWsR, roiNames, 'ad3ecR', 'vad16ecR', 'dlw', 1, 'ranksum');
%}

    % change Z score
%{
    cnDLWs = calcZScores(cnDLWs);
    adDLWs = calcZScores(adDLWs);
    vadDLWs = calcZScores(vadDLWs);
    vad2DLWs = calcZScores(vad2DLWs);
    vad3DLWs = calcZScores(vad3DLWs);
    vad4DLWs = calcZScores(vad4DLWs);
    vad5DLWs = calcZScores(vad5DLWs);
    vad6DLWs = calcZScores(vad6DLWs);
    cnDLWnss = calcZScores(cnDLWnss);
    adDLWnss = calcZScores(adDLWnss);
    vad6DLWnss = calcZScores(vad6DLWnss);
%}
    % plot correlation and cos similarity
    algNum = 25;
    meanCnDLW = nanmean(cnDLWs,3);
    meanAdDLW = nanmean(adDLWs,3);
    meanVadDLW = nanmean(vadDLWs,3);
    meanVad2DLW = nanmean(vad2DLWs,3);
%    meanVad3DLW = nanmean(vad3DLWs,3);
%    meanVad4DLW = nanmean(vad4DLWs,3);
%    meanVad5DLW = nanmean(vad5DLWs,3);
%    meanVad6DLW = nanmean(vad6DLWs,3);
    meanVad7DLW = nanmean(vad7DLWs,3);
    meanVad8DLW = nanmean(vad8DLWs,3);
    meanVad9DLW = nanmean(vad9DLWs,3);
    meanVad10DLW = nanmean(vad10DLWs,3);
    meanVad12DLW = nanmean(vad12DLWs,3);
    meanVad19DLW = nanmean(vad19DLWs,3);
%    figure; cnadDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanAdDLW);
%    figure; cnvadDLWr = plotTwoSignalsCorrelation(meanCnDLW, meanVadDLW);
    figure; advadDLWr = plotTwoSignalsCorrelation(meanAdDLW, meanVadDLW + nanx);
    figure; advadDLWr2 = plotTwoSignalsCorrelation(meanAdDLW, meanVad2DLW);
%    figure; advadDLWr3 = plotTwoSignalsCorrelation(meanAdDLW, meanVad3DLW);
%    figure; advadDLWr4 = plotTwoSignalsCorrelation(meanAdDLW, meanVad4DLW);
%    figure; advadDLWr5 = plotTwoSignalsCorrelation(meanAdDLW, meanVad5DLW);
%    figure; advadDLWr6 = plotTwoSignalsCorrelation(meanAdDLW, meanVad6DLW);
%    figure; advadDLWr7 = plotTwoSignalsCorrelation(meanAdDLW, meanVad7DLW + nanx);
%    figure; advadDLWr8 = plotTwoSignalsCorrelation(meanAdDLW, meanVad8DLW + nanx);
%    figure; advadDLWr9 = plotTwoSignalsCorrelation(meanAdDLW, meanVad9DLW + nanx);
%    figure; advadDLWr10 = plotTwoSignalsCorrelation(meanAdDLW, meanVad10DLW + nanx);
%    figure; advadDLWr12 = plotTwoSignalsCorrelation(meanAdDLW, meanVad12DLW + nanx);
    figure; advadDLWr19 = plotTwoSignalsCorrelation(meanAdDLW, meanVad19DLW + nanx);
    cosSim = zeros(algNum,1);
    cosSim(1) = getCosSimilarity(meanCnDLW, meanAdDLW);
    cosSim(2) = getCosSimilarity(meanCnDLW, meanVadDLW);
    cosSim(3) = getCosSimilarity(meanAdDLW, meanVadDLW);
    cosSim(4) = getCosSimilarity(meanCnDLW, meanVad2DLW);
    cosSim(5) = getCosSimilarity(meanAdDLW, meanVad2DLW);
%    cosSim(6) = getCosSimilarity(meanCnDLW, meanVad3DLW);
%    cosSim(7) = getCosSimilarity(meanAdDLW, meanVad3DLW);
%    cosSim(8) = getCosSimilarity(meanCnDLW, meanVad4DLW);
%    cosSim(9) = getCosSimilarity(meanAdDLW, meanVad4DLW);
%    cosSim(10) = getCosSimilarity(meanCnDLW, meanVad5DLW);
%    cosSim(11) = getCosSimilarity(meanAdDLW, meanVad5DLW);
%    cosSim(12) = getCosSimilarity(meanCnDLW, meanVad6DLW);
%    cosSim(13) = getCosSimilarity(meanAdDLW, meanVad6DLW);
%    cosSim(14) = getCosSimilarity(meanCnDLW, meanVad7DLW);
%    cosSim(15) = getCosSimilarity(meanAdDLW, meanVad7DLW);
%    cosSim(16) = getCosSimilarity(meanCnDLW, meanVad8DLW);
%    cosSim(17) = getCosSimilarity(meanAdDLW, meanVad8DLW);
    cosSim(18) = getCosSimilarity(meanCnDLW, meanVad9DLW);
    cosSim(19) = getCosSimilarity(meanAdDLW, meanVad9DLW);
    cosSim(20) = getCosSimilarity(meanCnDLW, meanVad10DLW);
    cosSim(21) = getCosSimilarity(meanAdDLW, meanVad10DLW);
    cosSim(22) = getCosSimilarity(meanCnDLW, meanVad12DLW);
    cosSim(23) = getCosSimilarity(meanAdDLW, meanVad12DLW);
    cosSim(24) = getCosSimilarity(meanCnDLW, meanVad19DLW);
    cosSim(25) = getCosSimilarity(meanAdDLW, meanVad19DLW);
    X = categorical({'cn-ad','cn-vad','ad-vad','cn-vad2','ad-vad2','cn-vad3','ad-vad3','cn-vad4','ad-vad4','cn-vad5','ad-vad5',...
        'cn-vad6','ad-vad6','cn-vad7','ad-vad7','cn-vad8','ad-vad8','cn-vad9','ad-vad9','cn-vad10','ad-vad10','cn-vad11','ad-vad11',...
        'cn-vad19','ad-vad19'});
    figure; bar(X, cosSim);
    title('cos similarity between CN and AD by each algorithm');

    % normality test
%{
    cnDLWsNt = calculateAlzNormalityTest(cnDLWs, roiNames, 'cnec', 'dlw');
    adDLWsNt = calculateAlzNormalityTest(adDLWs, roiNames, 'adec', 'dlw');
    vadDLWsNt = calculateAlzNormalityTest(vadDLWs, roiNames, 'vadec', 'dlw');
    vad7DLWsNt = calculateAlzNormalityTest(vad7DLWs, roiNames, 'vad7ec', 'dlw');
    cnnsDLWsNt = calculateAlzNormalityTest(cnDLWnss, roiNames, 'cnns', 'dlw');
    adnsDLWsNt = calculateAlzNormalityTest(adDLWnss, roiNames, 'adns', 'dlw');
    vadnsDLWsNt = calculateAlzNormalityTest(vadDLWnss, roiNames, 'vadns', 'dlw');
    vad3nsDLWsNt = calculateAlzNormalityTest(vad3DLWnss, roiNames, 'vad3ns', 'dlw');
    vad4nsDLWsNt = calculateAlzNormalityTest(vad4DLWnss, roiNames, 'vad4ns', 'dlw');
%}
    % compalizon test (Wilcoxon, Mann?Whitney U test)
    [cnadDLWsUt, cnadDLWsUtP, cnadDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, adDLWs, roiNames, 'cnec', 'adec', 'dlw');
    [cnvadDLWsUt, cnvadDLWsUtP, cnvadDLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vadDLWs, roiNames, 'cnec', 'vadec', 'dlw');
    [advadDLWsUt, advadDLWsUtP, advadDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vadDLWs, roiNames, 'adec', 'vadec', 'dlw');
%    [~, ~, ~] = calculateAlzWilcoxonTest(adDLWsR, vadDLWsR, roiNames, 'adecR', 'vadecR', 'dlw');
%    [cnvad2DLWsUt, cnvad2DLWsUtP, cnvad2DLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vad2DLWs, roiNames, 'cnec', 'vad2ec', 'dlw');
    [advad2DLWsUt, advad2DLWsUtP, advad2DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad2DLWs, roiNames, 'adec', 'vad2ec', 'dlw');
%    [~, ~, ~] = calculateAlzWilcoxonTest(adDLWsR, vad2DLWsR, roiNames, 'adecR', 'vad2ecR', 'dlw');
%    [cnvad7DLWsUt, cnvad7DLWsUtP, cnvad7DLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vad7DLWs, roiNames, 'cnec', 'vad7ec', 'dlw');
%    [advad7DLWsUt, advad7DLWsUtP, advad7DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad7DLWs, roiNames, 'adec', 'vad7ec', 'dlw');
%    [~, ~, ~] = calculateAlzWilcoxonTest(adDLWsR, vad7DLWsR, roiNames, 'adecR', 'vad7ecR', 'dlw');
%    [cnvad8DLWsUt, cnvad8DLWsUtP, cnvad8DLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vad8DLWs, roiNames, 'cnec', 'vad8ec', 'dlw');
%    [advad8DLWsUt, advad8DLWsUtP, advad8DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad8DLWs, roiNames, 'adec', 'vad8ec', 'dlw');
    [cnvad9DLWsUt, cnvad9DLWsUtP, cnvad9DLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vad9DLWs, roiNames, 'cnec', 'vad9ec', 'dlw');
    [advad9DLWsUt, advad9DLWsUtP, advad9DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad9DLWs, roiNames, 'adec', 'vad9ec', 'dlw');
%    [advad2bDLWsUt, advad2bDLWsUtP, advad2bDLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad2bDLWs, roiNames, 'adec', 'vad2bec', 'dlw');
%    [advad3DLWsUt, advad3DLWsUtP, advad3DLWsUtP2] = calculateAlzWilcoxonTest(cnDLWs, vad3DLWs, roiNames, 'cnec', 'vad3ec', 'dlw');
%    [advad3DLWsUt, advad3DLWsUtP, advad3DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad3DLWs, roiNames, 'adec', 'vad3ec', 'dlw');
%    [advad5DLWsUt, advad5DLWsUtP, advad5DLWsUtP2] = calculateAlzWilcoxonTest(adDLWs, vad5DLWs, roiNames, 'adec', 'vad5ec', 'dlw');
    [cnadDLWnssUt, cnadDLWnssUtP, cnadDLWnssUtP2] = calculateAlzWilcoxonTest(cnDLWnss, adDLWnss, roiNames, 'cnns', 'adns', 'dlw');
    [cnvadDLWnssUt, cnvadDLWnssUtP, cnvadDLWnssUtP2] = calculateAlzWilcoxonTest(cnDLWnss, vadDLWnss, roiNames, 'cnns', 'vadns', 'dlw');
    [advadDLWnssUt, advadDLWnssUtP, advadDLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vadDLWnss, roiNames, 'adns', 'vadns', 'dlw');
    [advad2DLWnssUt, advad2DLWnssUtP, advad2DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad2DLWnss, roiNames, 'adns', 'vad2ns', 'dlw');
%    [advad7DLWnssUt, advad7DLWnssUtP, advad7DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad7DLWnss, roiNames, 'adns', 'vad7ns', 'dlw');
    [advad9DLWnssUt, advad9DLWnssUtP, advad9DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad9DLWnss, roiNames, 'adns', 'vad9ns', 'dlw');
%    [advad3DLWnssUt, advad3DLWnssUtP, advad3DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad3DLWnss, roiNames, 'adns', 'vad3ns', 'dlw');
%    [advad4DLWnssUt, advad4DLWnssUtP, advad4DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad4DLWnss, roiNames, 'adns', 'vad4ns', 'dlw');
%    [advad5DLWnssUt, advad5DLWnssUtP, advad5DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad5DLWnss, roiNames, 'adns', 'vad5ns', 'dlw');
%    [advad6DLWnssUt, advad6DLWnssUtP, advad6DLWnssUtP2] = calculateAlzWilcoxonTest(adDLWnss, vad6DLWnss, roiNames, 'adns', 'vad6ns', 'dlw');
    
%{
    % using minimum 100 p-value relations. perform 5-fold cross validation.
    topNum = 100;
    sigTh = 2;
    N = 1;

    dlwAUC = zeros(1,N);
    dlwvAUC = zeros(1,N);
    dlwROC = cell(N,2);
    dlwvROC = cell(N,2);

    sigCntCN = cell(N,algNum);
    sigCntAD = cell(N,algNum);
    for k=1:N
        i = 1;
        % check sigma of healthy subject
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs, adDLWs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnadDLWsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlwROC{k,1}, dlwROC{k,2}, dlwAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC

        % check sigma of healthy subject
        i = i + 1;
        [control, target, meanTarget, stdTarget, meanControl] = getkFoldDataSet(cnDLWs, vadDLWs, k, N);         % replece cn*s, ad*s
        [B, I, X] = sortAndPairPValues(control, target, cnvadDLWsUtP, topNum);                                  % replace cnad*sUtP
        sigCntCN{k,i} = calcAlzSigmaSubjects(control, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        sigCntAD{k,i} = calcAlzSigmaSubjects(target, meanTarget, stdTarget, meanControl, I, topNum, sigTh);
        [dlwvROC{k,1}, dlwvROC{k,2}, dlwvAUC(k)] = calcAlzROCcurve(sigCntCN{k,i}, sigCntAD{k,i}, topNum);         % replace *ROC, *AUC
    end
    figure; boxplot(X);

    % save result
    fname = ['results/adsim-cn-ad-roi' num2str(132) '-result.mat'];
    save(fname, 'cosSim', 'dlwAUC','dlwvAUC','dlwROC','dlwvROC', 'sigCntCN', 'sigCntAD');
    mean(dlwAUC) % show result AUC
    mean(dlwvAUC) % show result AUC

    % plot ROC curve
    figure;
    hold on;
    plotErrorROCcurve(dlwROC, N, [0.2,0.2,0.2]); % TODO:
    plotErrorROCcurve(dlwvROC, N, [0.6,0.2,0.2]); % TODO:
    plotAverageROCcurve(dlwROC, N, '-', [0.2,0.2,0.2],1.0);
    plotAverageROCcurve(dlwvROC, N, '-', [0.6,0.2,0.2],1.0);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title(['averaged ROC curve']);
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
%}
end

% ==================================================================================================================

function [r1m, r2m, r3m, h1c, p1m] = retrainDLCMAndECmultiPattern(cnSignals, adDLWs, vad19DLWs, vad19DLWnss, vad19Zij, vad19DLWsR, cnS2, cnIS2, roiNames, group)
    ROINUM = size(cnS2,1);
    cnSbjNum = size(cnS2,3);
    nanx = eye(ROINUM);
    nanx(nanx==1) = NaN;

    R = 4; %ROINUM;
    JMAX = 7;
    k1 = floor(101/20)+1;
    r1 = zeros(JMAX+1,k1,R);
    r2 = zeros(JMAX+1,k1,R,cnSbjNum);
    r3 = zeros(JMAX+1,k1,cnSbjNum);
    h1 = zeros(JMAX+1,k1,ROINUM,ROINUM);
    p1 = zeros(JMAX+1,k1,ROINUM,ROINUM);
    h1c = zeros(JMAX+1,k1);
    for i=0:0
        for j=0:6
            for k=1:20:101
                vad20DLWnss = vad19DLWnss;
                if strcmp(group, 'vad24')
                    vad20DLWnss(:,1:ROINUM,:) = vad19Zij - vad19DLWsR * j * 0.2;
                    vad20DLWnss = [vad20DLWnss(:,1:ROINUM,:) repmat(vad20DLWnss(:,ROINUM+1,:),[1 k 1]) vad20DLWnss(:,ROINUM+2:end,:)];
                    cnS20 = [cnS2(:,1:ROINUM,:) repmat(cnS2(:,ROINUM+1,:),[1 k 1]) cnS2(:,ROINUM+2:end,:)];
                    cnIS20 = [cnIS2(:,1:ROINUM,:) repmat(cnIS2(:,ROINUM+1,:),[1 k 1]) cnIS2(:,ROINUM+2:end,:)];
                else
                    vad20DLWnss(:,2:ROINUM+1,:) = vad19Zij - vad19DLWsR * j * 0.2;
                    vad20DLWnss = [repmat(vad20DLWnss(:,1,:),[1 k 1]) vad20DLWnss(:,2:end,:)];
                    cnS20 = [repmat(cnS2(:,1,:),[1 k 1]) cnS2(:,2:end,:)];
                    cnIS20 = [repmat(cnIS2(:,1,:),[1 k 1]) cnIS2(:,2:end,:)];
                end

                vad21name = [group '-' num2str(i) '-' num2str(j) '-' num2str(k) 'ns'];
                vad21ecname = [group '-' num2str(i) '-' num2str(j) '-' num2str(k) 'ec'];
                [vad21DLWs, meanVad21DLWns, stdVad21DLWns] = retrainDLCMAndEC(vad20DLWnss, cnS20, cnIS20, roiNames, vad21name);
                [vad21bDLWs, vad21DLWnss] = calculateNodeSignals(cnSignals, cnS20, cnIS20, roiNames, vad21name, 'dlw');

                k1 = floor(k/20)+1;
                if strcmp(group, 'vad24')
                    vad19bDLWnss = [vad19DLWnss(:,1:ROINUM,:) repmat(vad19DLWnss(:,ROINUM+1,:),[1 k 1]) vad19DLWnss(:,ROINUM+2:end,:)];
                else
                    vad19bDLWnss = [repmat(vad19DLWnss(:,1,:),[1 k 1]) vad19DLWnss(:,2:end,:)];
                end
                for b=1:R
                    r1(j+1,k1,b) = corr2(squeeze(vad19bDLWnss(b,1,:)), squeeze(vad21DLWnss(b,1,:)));
%                    figure; hold on; plot([0.6 1.1], [0.6 1.1],':','Color',[0.5 0.5 0.5]); title(['nss corr: ' vad21name ' row=' num2str(b)]);
                    for a=1:cnSbjNum
%                        plotTwoSignalsCorrelation(vad19bDLWnss(b,1,a), vad21DLWnss(b,1,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.5], 'd', 8);
%                        plotTwoSignalsCorrelation(vad19bDLWnss(b,k+1:k+66,a), vad21DLWnss(b,k+1:k+66,a), [0.1*mod(a,10) 0.2*ceil(a/10) 0.8]);
                        r2(j+1,k1,b,a) = corr2(vad19bDLWnss(b,k+1:k+ROINUM,a), vad21DLWnss(b,k+1:k+ROINUM,a));
                    end; hold off;
                end
%                figure; hold on; plot([0 0.2], [0 0.2],':','Color',[0.5 0.5 0.5]); title(['ec corr: ' vad21ecname ' row=' num2str(b)]);
                for a=1:cnSbjNum
                    X = vad19DLWs(1:R,1:R,a)+nanx(1:R,1:R);
                    Y = vad21bDLWs(1:R,1:R,a);
%                    plotTwoSignalsCorrelation(X, Y, [0.1*mod(a,10) 0.2*ceil(a/10) 0.5]);
                    r3(j+1,k1,a) = corr2(X(~isnan(X(:))), Y(~isnan(Y(:))));
                end; hold off;
%                calculateAlzWilcoxonTest(vad19bDLWnss, vad21DLWnss, roiNames, 'vad19ns', vad21name, 'dlw', 1, 'ranksum');
                [h1(j+1,k1,:,:), p1(j+1,k1,:,:), ~] = calculateAlzWilcoxonTest(adDLWs, vad21bDLWs, roiNames, 'adec', vad21ecname, 'dlw', 1, 'ranksum', 1);
                h1c(j+1,k1) = length(find(h1(j+1,k1,1:R,:)>0));
            end
        end
    end
    r1m = nanmean(r1,3);
    r2m = nanmean(nanmean(r2,4),3);
    r3m = nanmean(r3,3);
    p1m = nanmean(nanmean(p1(:,:,1:R,:),4),3);
end

function sigEC = sigmaEC(EC)
    m = nanmean(EC(:));
    s = nanstd(EC(:),1);
    sigEC = (EC - m) / s;
end

function [weights, meanWeights, stdWeights] = retrainDLCMAndEC(teachSignals, nodeSignals, exSignals, roiNames, group)
    ROWNUM = size(teachSignals,1);
    COLNUM = size(teachSignals,2);
    sbjNum = size(teachSignals,3);
    weights = zeros(ROWNUM, ROWNUM, sbjNum);

    outfName = ['results/adsim-retrain-' group '-roi' num2str(ROWNUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        % init params
        sigLen = size(nodeSignals,2);
        si = nodeSignals;
        inSignal = exSignals;
        inControl = eye(ROWNUM);

        for i=1:sbjNum
            dlcmName = ['results/adsim-dlcm-' group '-roi' num2str(ROWNUM) '-net' num2str(i) '.mat'];
            if exist(dlcmName, 'file')
                load(dlcmName);
            else
                % init DLCM network
                netDLCM = initDlcmNetwork(si, inSignal, [], inControl);

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
                for j=1:4 %ROWNUM
                    nodeTeach = teachSignals(j,1:end,i);
                    nodeInput = [nodeSignals; exSignals];
                    if ~isempty(inControl)
                        filter = repmat(inControl(i,:).', 1, size(nodeInput,2));
                        nodeInput(ROWNUM+1:end,:) = nodeInput(ROWNUM+1:end,:) .* filter;
                    end
                    idx = find(isnan(nodeTeach));
                    nodeTeach(:,idx) = [];
                    nodeInput(:,idx) = [];
                    [netDLCM.nodeNetwork{j}, netDLCM.trainInfo{j}] = trainNetwork(nodeInput, nodeTeach, netDLCM.nodeLayers{j}, options);
                    disp(['virtual alzheimer (' group ') training node ' num2str(i) '-' num2str(j) ' rmse=' num2str(netDLCM.trainInfo{j}.TrainingRMSE(maxEpochs))]);
                end

                save(dlcmName, 'netDLCM', 'si', 'inSignal', 'inControl', 'options');
            end

            % recalculate EC
            weights(:,:,i) = calcDlcmEC(netDLCM, [], inControl);
        end
        save(outfName, 'weights', 'roiNames');
    end
    meanWeights = nanmean(weights, 3);
    stdWeights = nanstd(weights, 1, 3);
end

function vadSignals = calculateVirtualADSignals3(signals, roiNames, cnECs, adECs, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    vadSignals = {};
    
    outfName = ['results/adsim-signal-vad-roi' num2str(ROINUM) '-' algorithm '.mat'];
    if exist(outfName, 'file')
        load(outfName);
        return;
    end
    
    cnOutECs = nansum(nanmean(cnECs, 3), 1);
    adOutECs = nansum(nanmean(adECs, 3), 1);
    transRate = adOutECs ./ cnOutECs;

    for i=1:length(signals)
        sbjSignals = signals{i};
        for j=1:ROINUM
            sbjSignals(j,:) = sbjSignals(j,:) * transRate(1,j);
        end
        vadSignals{end+1} = sbjSignals;
    end
    save(outfName, 'vadSignals', 'roiNames');
end

function [vadSignals, vadDLWs, vadDLWnss] = calculateVirtualADSignals4(cnSignals, adSignals, roiNames, cnInSignals, cnInControls, group)
    ROWNUM = size(cnSignals{1},1);
    sigLen = size(cnSignals{1},2);
    cnNum = length(cnSignals);
    adNum = length(adSignals);
    outfName = ['results/adsim-all-' group '-roi' num2str(ROWNUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        allCnSignals = nan(ROWNUM, sigLen, cnNum);
        allVadSignals = nan(ROWNUM, sigLen, cnNum, adNum);
        sig = nan(cnNum);
        c = nan(cnNum);
        maxsi = nan(cnNum);
        minsi = nan(cnNum);
        for k=1:adNum
            disp(['generate virtual ad signals (type 4) : ' num2str(k)]);
            origName = ['results/ad-dlcm-ad-roi' num2str(ROWNUM) '-net' num2str(k) '.mat'];
            load(origName);
            for i=1:cnNum
                [allCnSignals(:,:,i), sig(i), c(i), maxsi(i), minsi(i)] = convert2SigmoidSignal(cnSignals{i});
                [Y, time] = predictDlcmNetwork(allCnSignals(:,:,i), cnInSignals(:,:,i), [], cnInControls(:,:,1), netDLCM);
                allVadSignals(:,:,i,k) = Y;
            end
        end
        save(outfName, 'allCnSignals', 'cnInSignals', 'sig', 'c', 'maxsi', 'minsi', 'allVadSignals');
    end
    % get mean of AD DLCM generated signals (type 4)
    meanVadSignals = nanmean(allVadSignals, 4);
    vadSignals = {};
    for i=1:cnNum
        %vadSignals{end+1} = convert2InvSigmoidSignal(meanVadSignals(:,:,i), sig(i), c(i), maxsi(i), minsi(i));
        vadSignals{end+1} = meanVadSignals(:,:,i);
        %plot(vadSignals{end});
    end
    [vadDLs, ~, ~] = calculateConnectivity(vadSignals, roiNames, group, 'dlcm', 1);
    [vadDLWs, ~, ~] = calculateConnectivity(vadSignals, roiNames, group, 'dlw', 1);
    [~, vadDLWnss, meanVadDLWns, stdVadDLWns, ~, ~] = calculateDistributions2(vadSignals, roiNames, group, 'dlw', group);
end

function [ECs, nodeSignals, meanSignals, stdSignals, inSignals, inControls] = calculateDistributions(signals, roiNames, group, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    outfName = ['results/adsim-' algorithm '-' group '-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        ECs = zeros(ROINUM, ROINUM, sbjNum);
        nodeSignals = zeros(ROINUM, ROINUM+1, sbjNum);
        inSignals = zeros(ROINUM, size(signals{1},2), sbjNum);
        inControls = zeros(ROINUM, ROINUM, sbjNum);
        for i=1:sbjNum
            switch(algorithm)
            case 'dlw'
                dlcmName = ['results/ad-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                load(dlcmName);
                [ec, ecNS] = calcDlcmEC(netDLCM, [], inControl);
            end
            ECs(:,:,i) = ec;
            nodeSignals(:,:,i) = ecNS;
            inSignals(:,:,i) = inSignal;
            inControls(:,:,i) = inControl;
        end
        save(outfName, 'ECs', 'nodeSignals', 'roiNames', 'inSignals', 'inControls');
    end
    meanSignals = nanmean(nodeSignals, 3);
    stdSignals = nanstd(nodeSignals, 1, 3);
    save(outfName, 'ECs', 'nodeSignals', 'meanSignals', 'stdSignals', 'roiNames', 'inSignals', 'inControls');
end

function [ECs, nodeSignals] = calculateNodeSignals(signals, S2, IS2, roiNames, group, algorithm)
    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    outfName = ['results/adsim-' algorithm '-' group '_ns-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        ECs = zeros(ROINUM, ROINUM, sbjNum);
        nodeSignals = zeros(ROINUM, size(S2, 2), sbjNum);
        for i=1:sbjNum
            switch(algorithm)
            case 'dlw'
                dlcmName = ['results/adsim-dlcm-' group '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                load(dlcmName);
                [Y, time] = predictDlcmNetwork(S2, IS2, [], inControl, netDLCM);
                ec = calcDlcmEC(netDLCM, [], inControl);
            end
            ECs(:,:,i) = ec;
            nodeSignals(:,:,i) = Y;
        end
        save(outfName, 'ECs', 'nodeSignals', 'roiNames', 'S2', 'IS2');
    end
end

function [ECs, nodeSignals, inSignals, inControls] = calculateNodeSignals2(signals, S2, IS2, roiNames, group, algorithm, prefix, orgGroup)
    % constant value
    ROINUM = size(signals{1},1);
    sbjNum = length(signals);

    outfName = ['results/adsim-' algorithm '-' group '_ns-roi' num2str(ROINUM) '.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        ECs = zeros(ROINUM, ROINUM, sbjNum);
        nodeSignals = zeros(ROINUM, size(S2, 2), sbjNum);
        if strcmp(prefix, 'adsim')
            inSignals = zeros(ROINUM, size(S2, 2), sbjNum);
        else
            inSignals = zeros(ROINUM, size(signals{1},2), sbjNum);
        end
        inControls = zeros(ROINUM, ROINUM, sbjNum);
        for i=1:sbjNum
            switch(algorithm)
            case 'dlw'
                dlcmName = ['results/' prefix '-dlcm-' orgGroup '-roi' num2str(ROINUM) '-net' num2str(i) '.mat'];
                load(dlcmName);
                [Y, time] = predictDlcmNetwork(S2, IS2, [], inControl, netDLCM);
                ec = calcDlcmEC(netDLCM, [], inControl);
            end
            ECs(:,:,i) = ec;
            nodeSignals(:,:,i) = Y;
            inSignals(:,:,i) = inSignal;
            inControls(:,:,i) = inControl;
        end
        save(outfName, 'ECs', 'nodeSignals', 'roiNames', 'inSignals', 'inControls', 'S2', 'IS2');
    end
end

function pat = repmatPat(repNum, ROINUM)
    patIn = [repmat(1,[repNum 1]); repmat(0,[repNum 1])];
    pat = repmat(patIn, [ceil(ROINUM/length(patIn)) 1]);
    pat = pat(1:ROINUM,:);
    pat = [pat, 1-pat, pat*0.5, (1-pat)*0.5];
end

function [ECs, nodeSignals, meanSignals, stdSignals, inSignals, inControls, S2, IS2] = calculateDistributions2(signals, roiNames, group, algorithm, orgGroup)
    % constant value
    ROINUM = size(signals{1},1);

    % generate node input signals pattern
    S2 = ones(ROINUM, ROINUM+1);
    S2(:,2:end) = S2(:, 2:end) - eye(ROINUM);
    S2 = [S2, zeros(ROINUM, 1), eye(ROINUM), eye(ROINUM)*0.5];
    IS2 = [ones(ROINUM, ROINUM+1), zeros(ROINUM, 1), zeros(ROINUM, ROINUM*2)];
    for i=1:7
        pat = repmatPat(2^(i-1), ROINUM);
        S2 = [S2, pat];
        IS2 = [IS2, pat];
    end

    % calculate node signals from input pattern
    [ECs, nodeSignals, inSignals, inControls] = calculateNodeSignals2(signals, S2, IS2, roiNames, group, algorithm, 'ad', orgGroup);

    meanSignals = nanmean(nodeSignals, 3);
    stdSignals = nanstd(nodeSignals, 1, 3);
end

function pat = repmatPat3(repNum, ROINUM)
    patIn = [repmat(1,[repNum 1]); repmat(0,[repNum 1])];
    pat = repmat(patIn, [ceil(ROINUM/length(patIn)) 1]);
    pat = pat(1:ROINUM,:);
    pat = [pat, 1-pat, pat*0.75, (1-pat)*0.75, pat*0.5, (1-pat)*0.5, pat*0.25, (1-pat)*0.25, ...
        1-pat*0.75, 1-(1-pat)*0.75, 1-pat*0.5, 1-(1-pat)*0.5, 1-pat*0.25, 1-(1-pat)*0.25];
end

function [ECs, nodeSignals, meanSignals, stdSignals, inSignals, inControls, S2, IS2] = calculateDistributions3(signals, roiNames, group, algorithm, orgGroup)
    % constant value
    ROINUM = size(signals{1},1);

    % generate node input signals pattern
    S2 = ones(ROINUM, ROINUM) - eye(ROINUM);
    S2 = [S2 ones(ROINUM, 1), ones(ROINUM, 1)*0.75, ones(ROINUM, 1)*0.5, ones(ROINUM, 1)*0.25, zeros(ROINUM, 1)];
    IS2 = [ones(ROINUM, ROINUM), ones(ROINUM, 1), ones(ROINUM, 1)*0.75, ones(ROINUM, 1)*0.5, ones(ROINUM, 1)*0.25, zeros(ROINUM, 1)];
    for i=1:7
        pat = repmatPat3(2^(i-1), ROINUM);
        S2 = [S2, pat];
        IS2 = [IS2, pat];
    end

    % calculate node signals from input pattern
    [ECs, nodeSignals, inSignals, inControls] = calculateNodeSignals2(signals, S2, IS2, roiNames, group, algorithm, 'ad', orgGroup);

    meanSignals = nanmean(nodeSignals, 3);
    stdSignals = nanstd(nodeSignals, 1, 3);
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

function [normalities, normalitiesP] = calculateAlzNormalityTest(ECs, roiNames, group, algorithm)
    % constant value
    ROWNUM = size(ECs,1);
    COLNUM = size(ECs,2);

    outfName = ['results/adsim-' algorithm '-' group '-roi' num2str(ROWNUM) '-normality.mat'];
    if exist(outfName, 'file')
        load(outfName);
    else
        normalities = nan(ROWNUM, COLNUM);
        normalitiesP = nan(ROWNUM, COLNUM);
        for i=1:ROWNUM
            for j=1:COLNUM
                if isnan(ECs(i,j,1)), continue; end
                x = squeeze(ECs(i,j,:));
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
