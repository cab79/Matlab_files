function S=eeglab_preprocess(S)
%% PREPROCESSING FOR CONTINUOUS EEGLAB .SET FILES

S.func = 'prep';



if 1
    % GET FILE LIST
    S.path.file = S.path.(S.func);
    S = getfilelist(S);

    % change to the input directory
    eval(sprintf('%s', ['cd(''' S.path.prep ''')']));

    % report if there are no such files
    if isempty(S.prep.filelist)
        error('No files found!\n');
    end
    
    % run though all files in a loop
    for f = S.prep.startfile:length(S.prep.filelist)
        filename = S.prep.filelist{f};

        fprintf('\nProcessing %s.\n\n', filename);

        % GET FILENAME PARTS
        [pth nme ext] = fileparts(filename); 
        fparts = strsplit(nme,'_');

        % LOAD DATA
        EEG = pop_loadset('filename',filename,'filepath',S.path.prep);

        %ADD CHANNEL LOCATIONS 
        if S.prep.chan.addloc && isempty(EEG.chanlocs(1).theta)
            EEG=pop_chanedit(EEG, 'lookup',S.path.locfile);
        end

        % SELECT TIME WINDOW TO ANALYSE
        if ~isempty(S.prep.cont.timewin)
            [~,~,i]=intersect(fparts,S.prep.select.conds);
            if ~isempty(S.prep.cont.timewin(i))
                EEG = pop_select(EEG,'time',S.prep.cont.timewin{i});
            end
        end

        % INTERPOLATE CHANNELS
        if ~isempty(S.prep.chan.interp) && S.prep.chan.interp~=0
            EEG= eeg_interp(EEG, S.prep.chan.interp, 'spherical');
        end

        % RE_REFERENCE CHANNELS
        if S.prep.chan.reref == 1
            % re-reference to the common average
            EEG = pop_reref( EEG, []); 
        elseif S.prep.chan.reref==2
            %% PrepPipeline - performs robust re-referencing
            % conventional re-referencing to the common average, prior to
            % cleaning could introduce artefact to otherwise clean channels.
            % Robust averaging as done here gets around this.
            params = struct();
            params.lineFrequencies = [];
            params.detrendType = 'high pass';
            params.detrendCutoff = 1;
            params.referenceType = 'robust';
            params.meanEstimateType = 'median';
            params.interpolationOrder = 'post-reference';
            params.keepFiltered = false;
            params.ignoreBoundaryEvents = true;
            [EEG, computationTimes] = prepPipeline(EEG, params); 
        end

        % NOTCH FILTER
        if ~isempty(S.prep.filter.notch) && iscell(S.prep.filter.notch)
            for i = 1:length(S.prep.filter.notch)
                EEG = pop_eegfiltnew(EEG,S.prep.filter.notch{i}(1),S.prep.filter.notch{i}(2),[],1); 
            end
        end

        % FILTER
        if ~isempty(S.prep.filter.incl)
            if S.prep.filter.incl(1)>0; EEG = pop_eegfiltnew( EEG, S.prep.filter.incl(1), 0, [], 0);end
            if S.prep.filter.incl(2)>0; EEG = pop_eegfiltnew( EEG, 0, S.prep.filter.incl(2), [], 0);end
        end

        % DOWNSAMPLE
        % do this after filtering if max freq to be analysed is more than half
        % of the downsample freq (Nyquist law)
        if S.prep.cont.downsample
            EEG = pop_resample(EEG, S.prep.cont.downsample);
        end

        % EPOCH
        %create epochs if markers exist
        try
            EEG = pop_epoch( EEG, S.prep.epoch.markers, S.prep.epoch.timewin);
        end

        % if markers don't exist, add markers and epoch
        if S.prep.epoch.addmarker && isempty(EEG.epoch)
            Sr = EEG.srate; % sampling rate of data
            Ndp = Sr*(S.prep.epoch.timewin(2)-S.prep.epoch.timewin(1)+1);% number of data points per epoch
            Tdp = size(EEG.data,2);% total number of data points in file
            Mep = floor(Tdp/Ndp);% max possible number of epochs
            for i = 1:Mep
                EEG.event(1,i).type = S.prep.epoch.markers{1};
                EEG.event(1,i).latency = Sr*S.prep.epoch.timewin(1)+(i-1)*Ndp+1;
                EEG.event(1,i).urevent = S.prep.epoch.markers{1};
            end
            EEG = pop_epoch( EEG, S.prep.epoch.markers(1), S.prep.epoch.timewin);
        end


        % LINEAR DETREND
        if S.prep.epoch.detrend
            for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
        end

        % REMOVE BASELINE
        if S.prep.epoch.rmbase
            EEG = pop_rmbase( EEG, [S.prep.epoch.basewin(1)*1000, S.prep.epoch.basewin(2)*1000]);
        end

        EEG = eeg_checkset( EEG );

        % SAVE
        sname_ext = 'epoched';
        sname = [nme '_' sname_ext '.' S.prep.fname.ext{:}];
        if ~exist(fullfile(S.path.prep,sname_ext),'dir')
            mkdir(fullfile(S.path.prep,sname_ext));
        end
        EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext)); 
    end

    % update last filenameparts
    for fnp = 1:length(S.prep.fname.parts)
        try
            for i = 1:length(S.prep.select.([S.prep.fname.parts{fnp} 's']))

                S.prep.select.([S.prep.fname.parts{fnp} 's']){i} = strrep(S.prep.select.([S.prep.fname.parts{fnp} 's']){i},'.','_');

                %S.prep.([S.prep.fname.parts{end} 's']){i} = [S.prep.([S.prep.fname.parts{end} 's']){i} '_' sname_ext];
            end
        end
    end
elseif exist(fullfile(S.path.prep,'epoched'),'dir')
    sname_ext = 'epoched';
end

%% combine data files
if 0%S.prep.combinefiles
    clear OUTEEG
    
    
    % GET FILE LIST
    S.path.file = fullfile(S.path.prep,sname_ext);
    S = getfilelist(S,sname_ext);
    
    loadpath = fullfile(S.path.prep,sname_ext);
    for s = 1:length(S.prep.select.subjects)
        if isempty(S.prep.select.sessions)
            S.prep.select.sessions = {''};
        end
        for a = 1:length(S.prep.select.sessions)
        
            % FIND THE FILES FOR THIS SUBJECT
            if ~isempty(S.prep.select.sessions{a})
                subfiles = S.prep.filelist(find(not(cellfun('isempty', strfind(S.prep.filelist,S.prep.select.subjects{s}))) & not(cellfun('isempty', strfind(S.prep.filelist,S.prep.select.sessions{a})))));
            else
                subfiles = S.prep.filelist(find(not(cellfun('isempty', strfind(S.prep.filelist,S.prep.select.subjects{s})))));
            end
            
            
            % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
            for f = 1:length(subfiles)

                % get filename parts
                file = subfiles{f};
                [pth nme ext] = fileparts(file); 
                fparts = strsplit(nme,'_');

                % create base name (without extension)
                basename='';
                for i = 1:length(fparts)
                    if i==1;c= '';else;c= '_';end
                    basename = [basename c fparts{i}];
                end

                % load and merge
                EEG = pop_loadset('filename',file,'filepath',loadpath);
                [EEG.epoch(:).file] = deal(basename);
                if exist('OUTEEG','var')
                    OUTEEG = pop_mergeset(OUTEEG, EEG);
                else
                    OUTEEG = EEG;
                    OUTEEG.fileinfo = [];
                end
                OUTEEG.fileinfo = [OUTEEG.fileinfo,{EEG.epoch.file}];
                clear EEG
            end
            EEG=OUTEEG;
            [EEG.epoch.file] = deal(EEG.fileinfo{:});
            clear OUTEEG
            EEG = eeg_checkset( EEG );

            % SAVE
            sname_ext = 'combined';
            if ~isempty(S.prep.select.sessions{a})
                S.prep.select.sessions{a} = ['allsessions_'];
            end
            sname = [S.prep.study{:} '_' S.prep.select.subjects{s} '_' S.prep.select.sessions{a} sname_ext '.' S.prep.fname.ext{:}];
            if ~exist(fullfile(S.path.prep,sname_ext),'dir')
                mkdir(fullfile(S.path.prep,sname_ext));
            end
            EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext));
        end
    end
elseif exist(fullfile(S.path.prep,'combined'),'dir')
    sname_ext = 'combined';
end

% FIND THE NEW DATA FILES
%files = dir(fullfile(S.path.prep,['*' S.runsubject '*' sname_ext]));
%files_ana = 1:length(files);

% NOISY TRIAL AND CHANNEL REJECTION USING FIELDTRIP
if 0%~isempty(S.prep.clean.FTrej) && iscell(S.prep.clean.FTrej)
    
    % GET FILE LIST
    S.path.file = fullfile(S.path.prep,sname_ext);
    if strcmp(sname_ext,'combined')
        S.prep.fname.parts = {'study','subject','session','suffix','ext'};
        S.prep.select.conds = {};
        S.prep.select.blocks = {};
    end
    S = getfilelist(S,sname_ext);
    
    loadpath = fullfile(S.path.prep,sname_ext);
    for f = S.prep.startfile:length(S.prep.filelist)
        file = S.prep.filelist{f};
        EEG = pop_loadset('filename',file,'filepath',loadpath);
        
        % REMOVE FLAT CHANNELS using 1/var
        if isfield(S.prep.clean,'flatchan') && S.prep.clean.flatchan.varthresh > 0
            S = flat_channel_reject(S,EEG);
            EEG = eeg_interp(EEG, S.prep.clean.flatchan.rejchan);
            if length(EEG.chanlocs)>EEG.nbchan
                EEG.chanlocs(S.prep.clean.flatchan.rejchan)=[];
            end
            EEG = pop_select(EEG, 'notrial', S.prep.clean.flatchan.rejtrial);
        end

        % strategy - only remove a very small number of very bad trials / chans
        % before ICA - do further cleaning after ICA
        for i = 1:length(S.prep.clean.FTrej)
            %try
                EEG = FTrejman(EEG,S.prep.clean.FTrej{i}); 
            %end
        end

        % SAVE
        [pth nme ext] = fileparts(file); 
        sname_ext = 'manrej';
        sname = [nme '_' sname_ext '.' S.prep.fname.ext{:}];
        if ~exist(fullfile(S.path.prep,sname_ext),'dir')
            mkdir(fullfile(S.path.prep,sname_ext));
        end
        EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext)); 
    end
elseif exist(fullfile(S.path.prep,'manrej'),'dir')
    sname_ext = 'manrej';
end

% ICA
if S.prep.clean.ICA
    % GET FILE LIST 
    S.path.file = fullfile(S.path.prep,sname_ext);
    S = getfilelist(S,'combined_manrej');
    
    loadpath = fullfile(S.path.prep,sname_ext);
    for f = S.prep.startfile:length(S.prep.filelist)
        file = S.prep.filelist{f};
        EEG = pop_loadset('filename',file,'filepath',loadpath);

        %RUN ICA 
        numcomp = numcompeig(EEG);
        EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',numcomp);

        % SAVE
        [pth nme ext] = fileparts(file); 
        sname_ext = 'ICA';
        sname = [nme '_' sname_ext '.' S.prep.fname.ext{:}];
        if ~exist(fullfile(S.path.prep,sname_ext),'dir')
            mkdir(fullfile(S.path.prep,sname_ext));
        end
        EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext));
    end
end