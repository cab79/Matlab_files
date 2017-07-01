function filtereeg(basename,low,high)

loadpaths

EEG = pop_loadset('filename',[basename '.set'],'filepath',filepath);

EEG = pop_eegfilt( EEG, low, high, [], [0]); % Highpass filter cutoff freq. 1Hz.
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', [basename 'filt']); % Save as new dataset.