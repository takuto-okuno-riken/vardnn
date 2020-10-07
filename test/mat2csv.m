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

