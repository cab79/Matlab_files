function S=eeglab_import(S)

S.func = 'import';

% GET FILE LIST
S.path.file = S.path.raw;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.raw ''')']));

% report if there are no such files
if isempty(S.(S.func).filelist)
    error('No files found to import!\n');
end

% indices of S.filelist for each subject
col_ind = find(ismember(S.(S.func).designmat(1,:),{'subject'}));
designtab = cell2table(S.(S.func).designmat(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
[~,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);


for sub = 1:length(uni_ind)
    
    % FIND THE FILES
    subfiles = S.(S.func).filelist(file_ind==uni_ind(sub));
    subdirs = S.(S.func).dirlist(file_ind==uni_ind(sub));

    % run though subject files in a loop and convert
    for f = 1:length(subfiles)
        filename = subfiles{f};
        dirname = subdirs{f};
        %if length(file) > 1
        %    error('Expected 1 recording file. Found %d.\n',length(file));
        %else
            %filename = file.name;
            fprintf('\nProcessing %s.\n\n', filename);
        %end

        switch S.(S.func).fname.ext{:}
            case 'vhdr' % Brainvision
                if isfield(S.import,'module') && any(strcmp(S.import.module,'fileio'))
                    EEG = pop_fileio(fullfile(dirname,filename));
                else
                    EEG = pop_loadbv_CAB(dirname,filename); % CAB version changes filesnames to be the same as those on disk
                end
            case 'cnt' % Neuroscan
                EEG = pop_loadcnt(fullfile(dirname,filename));
            case {'mff',''}
                % requires toolbox: mffimport2.0 (Srivas Chennu)
                EEG = pop_readegimff(fullfile(dirname,filename));
            case 'bdf'
                EEG = pop_biosig(fullfile(dirname,filename),[],'BDF'); % NOT TESTED. May need to add in channel range instead of []
            case 'csv'
                % NEEDS CODE
        end
        EEG = eeg_checkset(EEG);


        %add channel locations
        if S.(S.func).chan.addloc
            EEG=pop_chanedit(EEG, 'lookup',S.path.locfile);
        end


        % select chans
        S.(S.func).chanlocs=EEG.chanlocs;
        S=select_chans(S);
        fprintf('Removing excluded channels.\n');
        EEG = pop_select(EEG,'channel',S.(S.func).inclchan);
        S.(S.func).chanlocs = S.(S.func).chanlocs(S.(S.func).inclchan);

        [pth nme ext] = fileparts(filename); 
        savename = nme;
        EEG.setname = sprintf([S.(S.func).save.prefix{:} '%s' S.(S.func).save.suffix{:}],savename); % the output file is called: basename_orig
        EEG.filename = sprintf([S.(S.func).save.prefix{:} '%s' S.(S.func).save.suffix{:} '.set'],savename);
        EEG.filepath = fullfile(S.path.prep,'imported');
        
        if ~exist(EEG.filepath,'dir')
            mkdir(EEG.filepath)
        end

        fprintf('Saving %s%s.\n', EEG.filepath, EEG.filename);
        pop_saveset(EEG,'filename', EEG.filename, 'filepath', EEG.filepath);
    end
end

