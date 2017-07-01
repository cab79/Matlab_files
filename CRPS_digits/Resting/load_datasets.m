loadpaths
files = dir(fullfile(filepath,'P*100Hz.Exp3.set'));
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab 

for f = 1:length(files)
    filename = files(f).name;
    EEG = pop_loadset('filename',filename,'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); % copy it to ALLEEG
end

 eeglab redraw