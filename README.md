# Deep-Learning based Causal Model (DLCM) toolbox

## Introduction
DLCM is a powerful tool of data-driven analysis and simulation technique to estimate Effective Connectivity (EC).
Based on DLCM framework, two types of EC are defined, such as DLCM-EC and DLCM-GC to measure causal relation among multiple time-series data.
<div align="center">
<img src="data/figure1.jpg">
</div>

This toolbox includes several causal analysis algorithms, such as DLCM-EC, DLCM-GC, multivariate Granger Causality, pair-wised Granger Causality,
linear Transfer Entropy, Functional Connectivity (Correlation), Partial Correlation and Wavelet Coherence to estimate EC from multiple node signals.
<div align="center">
<img src="data/figure3b.jpg">
</div>

Command line tool could perform EC estimation with several causal analysis algorithms from node signals in csv file or mat file,
then show output causal relational matrix and save data in csv file.

## Requirements: software
* MATLAB R2019a or later
* Deep Learning Toolbox ver12.1 or later

## Command line tool
~~~
>> dlcm -h
usage: dlcm [options] filename.csv ...
  -e, --dlec          output DLCM Effective Connectivity matrix result (<filename>_dlec.csv)
  -d, --dlgc          output DLCM Granger Causality matrix result (<filename>_dlgc.csv)
  -m, --mvgc          output multivaliate Granger Causality matrix result (<filename>_mvgc.csv)
  -g, --pwgc          output pair-wised Granger Causality matrix result (<filename>_pwgc.csv)
  -t, --te            output (LINUE) Transfer Entropy matrix result (<filename>_te.csv)
  -f, --fc            output Functional Conectivity matrix result (<filename>_fc.csv)
  -p, --pc            output Partial Correlation matrix result (<filename>_pc.csv)
  -w, --wc            output Wavelet Coherence matrix result (<filename>_wc.csv)
  --pval              save P-value matrix of DLCM-GC, mvGC, pwGC, TE, FC and PC (<filename>_*_pval.csv)
  --fval alpha        save F-value with <alpha> matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_fval.csv, <filename>_*_fcrit.csv)
  --aic               save AIC matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_aic.csv)
  --bic               save BIC matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_bic.csv)
  --groundtruth files calculate ROC curve and save AUC of DLCM-EC, DLCM-GC, mvGC, pwGC, TE, FC, PC and WC (<filename>_*_auc.csv)
  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)
  --transopt num      signal transform option <num> (for type 1:centroid value)
  --format type       save file format <type> 0:csv, 1:mat(each), 2:mat(all) (default:0)
  --lag num           time lag <num> for mvGC, pwGC and TE (default:3)
  --ex files          DLCM exogenouse input signal <files> (file1.csv[:file2.csv:...])
  --nctrl files       DLCM node status control <files> (file1.csv[:file2.csv:...])
  --ectrl files       DLCM exogenous input control <files> (file1.csv[:file2.csv:...])
  --epoch num         DLCM training epoch number <num> (default:1000)
  --l2 num            DLCM training L2Regularization <num> (default:0.05)
  --roiname files     ROI names <files> (file1.csv[:file2.csv:...])
  --showsig           show node status signals of <filename>.csv
  --showex            show exogenous input signals of <file1>.csv
  --showmat           show result matrix of DLCM-EC, DLCM-GC, mvGC, pwGC, TE, FC, PC and WC
  --showcg            show circle graph of DLCM-EC, DLCM-GC, mvGC, pwGC, TE, FC, PC and WC
  --showroc           show ROC curve (by GroundTruth) of DLCM-EC, DLCM-GC, mvGC, pwGC, TE, FC, PC and WC
  -v, --version       show version number
  -h, --help          show command line help
~~~

## Example Results
<div align="center">
<img src="data/figure9b.jpg">
</div>

## Citing DLCM toolbox
If you find DLCM useful in your research, please consider citing:  
Takuto Okuno, Alexander Woodward,
["DLCM: A Data-Driven Deep-Learning Based Effective Connectivity Estimation Toolbox"](https://yahoo.com/), under reviewing.

