S.func = 'multitest';
close all
pth='C:\Data\CORE\eeg\ana\prep\ana';

files = dir(fullfile(pth,'*.set'));

select_times = [100,600];

for f = 1%:length(files)
    EEG = pop_loadset('filename',files(f).name,'filepath',pth);
    dp = dsearchn(EEG.times',select_times');
    S=MultiOutliers(S,EEG.data(:,dp(1):dp(2),:));
end