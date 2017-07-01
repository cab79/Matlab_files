function dataimport_pib(basename)

eval(sprintf(['cd \C:\EEGdata\']));


filepath=['C:\EEGdata\']; % scripts path


filenames = dir(sprintf('%s%s*', filepath, basename));

if isempty(filenames)
    error('No files found to import!\n');
end

mfffiles = filenames(logical(cell2mat({filenames.isdir})));
if length(mfffiles) > 1
    error('Expected 1 MFF recording file. Found %d.\n',length(mfffiles));
else
    filename = mfffiles.name;
    fprintf('\nProcessing %s.\n\n', filename);
    EEG = pop_readegimff(sprintf('%s%s', filepath, filename),'datatype','PIB');
end

EEG = eeg_checkset(EEG);

% lpfreq = 100;
% fprintf('Low-pass filtering below %dHz...\n',lpfreq);
% EEG = pop_eegfilt(EEG, 0, lpfreq, [], [0], 0, 0, 'fir1', 0);
% hpfreq = 1;
% fprintf('High-pass filtering above %dHz...\n',hpfreq);
% EEG = pop_eegfilt(EEG, hpfreq, 0, [], [0], 0, 0, 'fir1', 0);


EEG.setname = sprintf('%s_pib_orig',basename);
EEG.filename = sprintf('%s_pib_orig.set',basename);
EEG.filepath = filepath;

fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);

