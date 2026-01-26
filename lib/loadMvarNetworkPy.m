%%
% load mVAR network with python pickle
% (need to setup python environment, https://jp.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html
% input:
%  path            net data file path (saved by vneumodpy)

function net = loadMvarNetworkPy(path)
    pickle = py.importlib.import_module('pickle');
    fh = py.open([path '/list.dat'], 'rb');
    dat = pickle.load(fh);
    fh.close();
    net.nodeNum = double(dat{1});
    net.exNum = double(dat{3});
    net.sigLen = double(dat{2});
    net.lags = double(dat{5});

    fh = py.open([path '/residuals.dat'], 'rb');
    dat = pickle.load(fh);
    fh.close();
    r = cell(length(dat),1);
    for i = 1:length(dat)
        d = dat{i};  % we need this line. then cast.
        ri = single(d);
        r{i} = ri';
    end
    net.rvec = r;

    fh = py.open([path '/regress.dat'], 'rb');
    dat = pickle.load(fh);
    fh.close();
    b = cell(length(dat),1);
    for i = 1:length(dat)
        d = dat{i};  % we need this line. then cast.
        bi = double(d);
        b{i} = bi';
    end
    net.bvec = b;
end