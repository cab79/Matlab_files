function S=eeglab_import(S)

% GET FILE LIST
S.filepath = S.rawpath;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.rawpath ''')']));

% report if there are no such files
if isempty(S.filelist)
    error('No files found to import!\n');
end

% run though all files in a loop and convert
for f = 1:length(S.filelist)
    filename = S.filelist{f};
    %if length(file) > 1
    %    error('Expected 1 recording file. Found %d.\n',length(file));
    %else
        %filename = file.name;
        fprintf('\nProcessing %s.\n\n', filename);
    %end
    
    switch S.inputfileext
        case 'vhdr' % Brainvision
            EEG = pop_loadbv_CAB(S.rawpath,filename); % CAB version changes filesnames to be the same as those on disk
        case 'cnt' % Neuroscan
            EEG = pop_loadcnt(fullfile(S.rawpath,filename));
    end
    EEG = eeg_checkset(EEG);

    
    %add channel locations
    if ~isempty(S.chan.addloc)
        EEG=pop_chanedit(EEG, 'lookup',S.chan.addloc);
    end
    
    fprintf('Removing excluded channels.\n');
    EEG = pop_select(EEG,'nochannel',S.chan.excl);

    [pth nme ext] = fileparts(file.name); 
    savename = nme;
    EEG.setname = sprintf([S.saveprefix '%s' S.savesuffix],savename); % the output file is called: basename_orig
    EEG.filename = sprintf([S.saveprefix '%s' S.savesuffix '.set'],savename);
    EEG.filepath = S.setpath;
    
    fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
    pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
end

