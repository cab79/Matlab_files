function dataimport(filepath, basename)% use the very specific file name, and the whole of it!

% loadpaths % this is just another .m file with one line of code specifying the directory you want to be in. you can comment on it and just include something like: eval(sprintf(['cd /imaging/sg04/stanimira_p4/eeg_emg/p33_s2/...']));

eval(sprintf(['cd ' filepath]));

chanexcl = [1,8,14,17,21,25,32,38,43,44,48,49,56,63,64,68,69,73,74,81,82,88,89,94,95,99,107,113,114,119,120,121,125,126,127,128];

filenames = dir(sprintf('%s%s*', filepath, basename));

if isempty(filenames)
    error('No files found to import!\n');
end

for f = 1:length(filenames)
    mfffiles = filenames(f);
    if length(mfffiles) > 1
        error('Expected 1 MFF recording file. Found %d.\n',length(mfffiles));
    else
        filename = mfffiles.name;
        fprintf('\nProcessing %s.\n\n', filename);
        EEG = pop_readegimff(sprintf('%s%s', filepath, filename));
    end

    EEG = eeg_checkset(EEG);

    %%% preprocessing

    % samprate = 250; %downsaple if srate is higher that 250!
    % fprintf('Downsampling to %dHz.\n',samprate);
    % EEG = pop_resample(EEG,samprate);

    fprintf('Removing excluded channels.\n');
    EEG = pop_select(EEG,'nochannel',chanexcl);

    % lpfreq = 40;
    % fprintf('Low-pass filtering below %dHz...\n',lpfreq);
    % EEG = pop_eegfilt(EEG, 0, lpfreq, [], [0], 0, 0, 'fir1', 0);
    % hpfreq = 0.1;
    % fprintf('High-pass filtering above %dHz...\n',hpfreq);
    % EEG = pop_eegfilt(EEG, hpfreq, 0, [], [0], 0, 0, 'fir1', 0);

    EEG.setname = sprintf('%s_orig',basename); % the output file is called: basename_orig
    EEG.filename = sprintf('%s_orig.set',basename);
    EEG.filepath = filepath;

    fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
    pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
end
