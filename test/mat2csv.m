%% for DCM generated signal
%{
prefix = 'results/net-pat4-8x300-idx';
N = 8;
idxs = [1,2,6,7,8];

for i=1:length(idxs)
    for j=1:N
        fname = [prefix num2str(idxs(i)) '-' num2str(j) '.mat'];
        if ~exist(fname, 'file')
            continue;
        end
        % load
        load(fname);
        
        mat = data.';

        % output csv file
        outfname = [prefix num2str(idxs(i)) '-' num2str(j) '.csv'];
        T = array2table(mat);
        writetable(T,outfname,'WriteVariableNames',false);
        disp(['output csv file : ' outfname]);
    end
end
%}
%% for time comparison with DCM generated signal
prefix = 'results/net-timeD-';
N = 12;
idxs = [8,12,16];
reps = [20,10,10];

for i=1:length(idxs)
    n = idxs(i);
    m = reps(i);
    for k=1:m
        fname = [prefix num2str(n) '-' num2str(N) 'x' num2str(k) '.mat'];

        if ~exist(fname, 'file')
            continue;
        end
        % load
        load(fname);
        
        mat = data.';

        % output csv file
        outfname = [prefix num2str(n) '-' num2str(N) 'x' num2str(k) '.csv'];
        T = array2table(mat);
        writetable(T,outfname,'WriteVariableNames',false);
        disp(['output csv file : ' outfname]);
    end
end

%% for reduced wongwang generated signal
%{
N = 8;
idxs = [16,32,48,64,80,98];

for i=1:length(idxs)
    for j=1:N
        fname = ['data/tvb-wongwang' num2str(idxs(i)) 'x55scan-pat' num2str(i) '-64hz-' num2str(j) '.mat'];
        if ~exist(fname, 'file')
            continue;
        end
        % load
        load(fname);
        
        mat = data.';

        % output csv file
        outfname = ['data/tvb-wongwang' num2str(idxs(i)) 'x55scan-pat' num2str(i) '-64hz-' num2str(j) '.csv'];
        T = array2table(mat);
        writetable(T,outfname,'WriteVariableNames',false);
        disp(['output csv file : ' outfname]);
    end
end
%}

%% for ADNI2 fMRI data
%{
roiNum = 132;
group = {'cn', 'ad', 'mci'};

for i=1:length(group)
    fname = ['data/ad-signal-' group{i} '-roi' num2str(roiNum) '.mat'];
    load(fname);

    for j=1:length(signals)
        % output csv file
        outfname = ['data/ad-signal-' group{i} '-roi' num2str(roiNum) '-' num2str(j) '.csv'];
        mat = signals{j};
        T = array2table(mat.');
        writetable(T,outfname,'WriteVariableNames',false);
        disp(['output csv file : ' outfname]);
    end
end
%}