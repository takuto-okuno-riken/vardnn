%%
% Caluclate Wavelet Coherence
% returns mean Wavelet Cross Spectrum matrix (mWCS), full set of Wavelet
% Coherence (WCOH) and Cross Spectrums (WCS)
% input:
%  X     multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  isFullNode   return both node & exogenous causality matrix (optional)

function [mWCS, WCOH, WCS] = calcWaveletCoherence(X, exSignal, nodeControl, exControl, isFullNode)
    if nargin < 5
        isFullNode = 0;
    end
    if nargin < 4
        exControl = [];
    end
    if nargin < 3
        nodeControl = [];
    end
    if nargin < 2
        exSignal = [];
    end
    nodeNum = size(X,1);
    nodeInNum = nodeNum + size(exSignal,1);
    if isFullNode==0, nodeMax = nodeNum; else nodeMax = nodeInNum; end
    
    % set node input
    if ~isempty(exSignal)
        X = [X; exSignal];
    end
    
    mWCS = ones(nodeNum,nodeMax);
    WCOH = cell(nodeNum,nodeMax);
    WCS = cell(nodeNum,nodeMax);
    for i=1:nodeNum
        nodeInput = X;
        if ~isempty(nodeControl)
            filter = nodeControl(i,:).';
            nodeInput(1:nodeNum,:) = nodeInput(1:nodeNum,:) .* filter;
        end
        if ~isempty(exControl)
            filter = exControl(i,:).';
            nodeInput(nodeNum+1:end,:) = nodeInput(nodeNum+1:end,:) .* filter;
        end
        for j=1:nodeMax
            if i==j, continue; end
            [wcoh, wcs] = wcoherence(nodeInput(i,:), nodeInput(j,:));
            WCOH{i,j} = wcoh;
            WCS{i,j} = wcs;
            mWCS(i,j) = nanmean(real(wcs),'all');
        end
    end
end
