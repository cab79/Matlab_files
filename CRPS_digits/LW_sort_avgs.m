clear all
path = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\LW\Averages';
conds = {
    {'avg P','avg H'};
    {'left','right'};
    {'D1','M3','D5'};
    };
types = [2 1 1 2 1 1 1 2 1 2 1 1 1];
type_name = {'Acc','UnAcc'};
type_num = 2;

for n = 1:length(conds{1});
    if ~exist(fullfile(path,conds{1}{n}),'dir')
        mkdir(fullfile(path,conds{1}{n}));
    end
    movefile(fullfile(path,['*' conds{1}{n} '*']),fullfile(path,conds{1}{n}));
    for n2 = 1:length(conds{2});
        if ~exist(fullfile(path,conds{1}{n},conds{2}{n2}),'dir')
            mkdir(fullfile(path,conds{1}{n},conds{2}{n2}));
        end
        movefile(fullfile(path,conds{1}{n},['*' conds{2}{n2} '*']),fullfile(path,conds{1}{n},conds{2}{n2}));
        for n3 = 1:length(conds{3});
            if ~exist(fullfile(path,conds{1}{n},conds{2}{n2},conds{3}{n3}),'dir')
                mkdir(fullfile(path,conds{1}{n},conds{2}{n2},conds{3}{n3}));
            end
            movefile(fullfile(path,conds{1}{n},conds{2}{n2},['*' conds{3}{n3} '*']),fullfile(path,conds{1}{n},conds{2}{n2},conds{3}{n3}));
        end
    end
end
