% set global random stream and shuffle it
myStream=RandStream('mt19937ar');
RandStream.setGlobalStream(myStream);
rng('shuffle');

% DEM Structure: create random inputs
% -------------------------------------------------------------------------
T  = 1000;                            % number of observations (scans)
TR = 2;                               % repetition time or timing
n  = 30;                              % number of regions or nodes
t  = (1:T)*TR;                        % observation times
k  = 6;

% priors
% -------------------------------------------------------------------------
options.maxnodes   = 4;  % effective number of nodes, 4 is better than n

options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 1;
options.induced    = 1;

A   = ones(n,n);
B   = zeros(n,n,0);
C   = zeros(n,n);
D   = zeros(n,n,0);
pP  = spm_dcm_fmri_priors(A,B,C,D,options);

pP.C = eye(n,n);
pP.transit = randn(n,1)/16;

% simulate response to endogenous fluctuations
%==========================================================================

% integrate states
% -------------------------------------------------------------------------
U.dt = TR;
M.f  = 'spm_fx_fmri';
M.x  = sparse(n,5);
M.g   = 'spm_gx_fmri';

pP.A = rand(n,n)/4 - 0.125;

Uus = {};
si = [];
uu = [];
for i=1:k
    U.u = spm_rand_mar(T+50,n,1/2)/8;       % endogenous fluctuations
    Uus{end+1} = U.u;
    y    = spm_int_J(pP,M,U);            % integrate with observer
    y2   = y(51:end,:);
    u2   = U.u(51:end,:);
    si = [si; y2.'];
    uu = [uu; u2.'];
    figure; plot(U.u);
    figure; plot(y2);
end

save(['test/testTrain-rand' num2str(n) '-dcm.mat'], 'si', 'uu', 'pP', 'Uus', 'M', 'T', 'TR','n', 't');
