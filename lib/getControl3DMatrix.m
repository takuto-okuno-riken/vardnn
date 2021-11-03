%%
% get control 3D matrix

function [nodeControl,exControl,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags)
    if isempty(nodeControl) && size(nodeControl,1)==0 && size(nodeControl,2)==0
        nodeControl = ones(nodeNum,nodeNum);
    end
    if isempty(exControl) && size(exControl,1)==0 && size(exControl,2)==0
        exControl = ones(nodeNum,exNum);
    end
    if lags > 1
        if ndims(nodeControl) <= 2
            nodeControl = repmat(nodeControl,[1,1,lags]);
        end
        if ndims(exControl) <= 2
            exControl = repmat(exControl,[1,1,lags]);
        end
    end
    control = [nodeControl, exControl];
end
