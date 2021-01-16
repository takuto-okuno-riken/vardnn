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

Command line tool could perform EC estimation with several causal analysis algorithms from node signals in csv or mat file,
then show output causal relational matrix and save data in csv or mat file.

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
  --outpath           output files path (default:"results")
  --pval              save P-value matrix of DLCM-GC, mvGC, pwGC, TE, FC and PC (<filename>_*_pval.csv)
  --fval alpha        save F-value with <alpha> matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_fval.csv, <filename>_*_fcrit.csv)
  --aic               save AIC matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_aic.csv)
  --bic               save BIC matrix of DLCM-GC, mvGC, pwGC and TE (<filename>_*_bic.csv)
  --format type       save file format <type> 0:csv, 1:mat(each), 2:mat(all) (default:0)
  --groundtruth files calculate ROC curve and save AUC of DLCM-EC, DLCM-GC, mvGC, pwGC, TE, FC, PC and WC (<filename>_*_auc.csv)
  --transform type    input signal transform <type> 0:raw, 1:sigmoid (default:0)
  --transopt num      signal transform option <num> (for type 1:centroid value)
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

## Command line tool Demo
This demo inputs 8 nodes random signal and outputs FC, mvGC, DLCM-GC and DLCM-EC results csv files and matrix graphs.
(Copy and paste this command line. Demo data is included in DLCM toolbox.)
~~~
>> dlcm -e -d -f -m --showsig --showmat --transform 1 --epoch 100 data/signal8.csv
start training
start training whole DLCM network
training node 1
training node 2
training node 3
training node 4
training node 5
training node 6
training node 7
training node 8
finish training whole DLCM network! t = 4.3548s
DLCM training result : rsme=0.16263
output csv file : results/signal8_dlec.csv
output csv file : results/signal8_dlgc.csv
output csv file : results/signal8_mvgc.csv
output csv file : results/signal8_fc.csv
~~~
These are output graphs of dlcm command.
<div align="center">
<img src="data/rdmfig1.jpg">
</div>

___
DLCM can take exogenous input signals with control matrix.
~~~
>> dlcm -e -d --showmat --epoch 100 --transform 1 --ex data/signal8ex.csv --ectrl data/ctrleye.csv data/signal8.csv
...
output csv file : results/signal8_dlec.csv
output csv file : results/signal8_dlgc.csv
~~~
___
This demo inputs 32 nodes synthetic fMRI BOLD signals of .mat file and outputs FC, PC, mvGC, TE, DLCM-GC and DLCM-EC results.
Result matrices of EC, P-value, F-value, AIC and BIC are saved in ww32-1_&lt;algorithm&gt;_all.mat file.
~~~
>> dlcm -e -d -f -p -m -t --transform 1 --pval --lag 5 --epoch 500 --l2 0.1 --fval 0.05 --aic --bic --format 2 --roiname roi32.csv --showsig --showmat --showcg --showroc data/ww32-1.mat data/ww32-2.mat data/ww32-3.mat data/ww32-4.mat
start training
start training whole DLCM network
training node 1
training node 2
...
training node 31
training node 32
finish training whole DLCM network! t = 61.5208s
DLCM training result : rsme=0.017795
~~~
.mat file includes input data matrices.
| name | matrix | description |
|:---|:---|:---|
|X |&lt;nodes&gt; x &lt;length&gt;(double)|node signals|
|exSignal|&lt;exogenous nodes&gt; x &lt;length&gt;(double)|exogenous signals|
|nodeControl|&lt;nodes&gt; x &lt;nodes&gt;(logical)|node connection control matrix|
|exControl|&lt;nodes&gt; x &lt;nodes&gt;(logical)|exogenous node connection control matrix|
|groundTruth|&lt;nodes&gt; x &lt;nodes&gt;(logical)|ground truth of network connection for ROC curve|

Several graphs (node signals, result matrix, circle graph, ROC curve) of each algorithm are shown by dlcm command.
<div align="center">
<img src="data/rdmfig2.jpg">
</div>

## Example Results
Example results of causal relation matrix graphs of human fMRI signals (132 ROI).
(Generating brain connectome image is not included in dlcm command)
<div align="center">
<img src="data/figure9b.jpg">
</div>

## Citing DLCM toolbox
If you find DLCM useful in your research, please consider citing:  
Takuto Okuno, Alexander Woodward,
["DLCM: A Data-Driven Deep-Learning Based Effective Connectivity Estimation Toolbox"](https://yahoo.com/), under reviewing.

