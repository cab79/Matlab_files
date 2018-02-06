function import_eeglab(S)
% Requires input S, a structure with the following fields:
%S.loadpath = ''; %where the raw input data is stored
%S.savepath = ''; %where the .set output files will go
%S.inputfileext = ''; % file extension of input data
%S.outputfileprefix = ''; % suffix to add to output file, if needed
%S.outputfilesuffix = ''; % suffix to add to output file, if needed
%S.chanexcl = []; % exclude channels

% change to the input directory
eval(sprintf('%s', ['cd(''' S.loadpath ''')']));

% find files with this extension in the input directory
filenames = dir(fullfile(S.loadpath,['*' S.inputfileext]));

% report if there are no such files
if isempty(filenames)
    error('No files found to import!\n');
end

% run though all files in a loop and convert
for f = 1:length(filenames)
    file = filenames(f);
    if length(file) > 1
        error('Expected 1 recording file. Found %d.\n',length(file));
    else
        filename = file.name;
        fprintf('\nProcessing %s.\n\n', filename);
    end
    
    switch S.inputfileext
        case 'vhdr' % Brainvision
            EEG = pop_loadbv_CAB(S.loadpath,filename); % CAB version changes filesnames to be the same as those on disk
        case 'cnt' % Neuroscan
            EEG = pop_loadcnt(fullfile(S.loadpath,filename));
    end
    EEG = eeg_checkset(EEG);

    
    %add channel locations
    if ~isempty(S.addchanloc)
        EEG=pop_chanedit(EEG, 'lookup',S.addchanloc);
    end
    
    fprintf('Removing excluded channels.\n');
    EEG = pop_select(EEG,'nochannel',S.chanexcl);

    [pth nme ext] = fileparts(file.name); 
    savename = nme;
    EEG.setname = sprintf([S.outputfileprefix '%s' S.outputfilesuffix],savename); % the output file is called: basename_orig
    EEG.filename = sprintf([S.outputfileprefix '%s' S.outputfilesuffix '.set'],savename);
    EEG.filepath = S.savepath;
    
    fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
    pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
end

