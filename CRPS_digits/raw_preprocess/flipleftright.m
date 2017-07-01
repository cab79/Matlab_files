EEG = flipchan(EEG);
EEG.saved = 'no';
[ALLEEG EEG index] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
fprintf('Channel laterality flipped.\n');