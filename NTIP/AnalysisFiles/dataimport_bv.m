function dataimport_bv(filepath, basename, parts, combine, savepath)% use the very specific file name, and the whole of it!

eval(sprintf('%s', ['cd(''' filepath ''')']));

chanexcl = [31,32];

filenames = dir(fullfile(filepath, basename));

if isempty(filenames)
    error('No files found to import!\n');
end

basesuff = '';
for f = 1:length(parts)
    %if ~strcmp(filenames(f).name(length(basename)+1),num2str(parts(f)));
    %    continue
    %end
    files = filenames(f);
    if length(files) > 1
        error('Expected 1 recording file. Found %d.\n',length(mfffiles));
    else
        filename = files.name;
        fprintf('\nProcessing %s.\n\n', filename);
        EEG = pop_loadbv(filepath,filename);
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

    [pth nme ext] = fileparts(files.name); 
    if length(parts)>1
        savename = [nme basesuff '_' num2str(parts(f))];
    else
        savename = nme;
    end
    
    EEG.setname = sprintf('%s_orig',savename); % the output file is called: basename_orig
    EEG.filename = sprintf('%s_orig.set',savename);
    EEG.filepath = savepath;
    
    if combine
        if exist('OUTEEG','var')
            OUTEEG = pop_mergeset(OUTEEG, EEG);
            if f==idx(end)
                EEG=OUTEEG;
                fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
                pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
            end
        else
            OUTEEG = EEG;
        end
        clear EEG
    else
        fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
        pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
    end
end

