epochlength = 30; %sec

events = (0:epochlength:EEG.xmax)';
events = cat(2,repmat({'EVNT'},length(events),1),num2cell(events));

EEG = pop_importevent(EEG,'event',events,'fields',{'type','latency'});
EEG = eeg_checkset(EEG,'makeur');
EEG = eeg_checkset(EEG,'eventconsistency');

fprintf('\nSegmenting into %d sec epochs.\n',epochlength);
epochEEG = pop_epoch(EEG,{'EVNT'},[0 epochlength]);

epochEEG = pop_rmbase(epochEEG,[]);

epochEEG = eeg_checkset(epochEEG);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, epochEEG);

eeglab redraw