function S=eeglab_preprocess(S)
%% PREPROCESSING FOR CONTINUOUS EEGLAB .SET FILES

% GET FILE LIST
S.filepath = S.setpath;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.setpath ''')']));

% report if there are no such files
if isempty(S.filelist)
    error('No files found!\n');
end

% run though all files in a loop
for f = 1:length(S.filelist)
    filename = S.filelist{f};
    
    fprintf('\nProcessing %s.\n\n', filename);

    % GET FILENAME PARTS
    [pth nme ext] = fileparts(filename); 
    fparts = strsplit(nme,'_');

    % LOAD DATA
    EEG = pop_loadset('filename',filename,'filepath',S.setpath);

    %ADD CHANNEL LOCATIONS 
    if ~isempty(S.chan.addloc) && ischar(S.chan.addloc) && isempty(EEG.chanlocs(1).theta)
        EEG=pop_chanedit(EEG, 'lookup',S.chan.addloc);
    end

    % SELECT TIME WINDOW TO ANALYSE
    if ~isempty(S.cont.timewin)
        [~,~,i]=intersect(fparts,S.conds);
        if ~isempty(S.cont.timewin(i))
            EEG = pop_select(EEG,'time',S.cont.timewin{i});
        end
    end

    % INTERPOLATE CHANNELS
    if ~isempty(S.chan.interp) && S.chan.interp~=0
        EEG= eeg_interp(EEG, S.chan.interp, 'spherical');
    end

    % RE_REFERENCE CHANNELS
    if S.chan.reref == 1
        % re-reference to the common average
        EEG = pop_reref( EEG, []); 
    elseif S.chan.reref==2
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
    if ~isempty(S.filter.notch) && iscell(S.filter.notch)
        for i = 1:length(S.filter.notch)
            EEG = pop_eegfiltnew(EEG,S.filter.notch{i}(1),S.filter.notch{i}(2),[],1); 
        end
    end

    % FILTER
    if S.filter.incl(1)>0; EEG = pop_eegfiltnew( EEG, S.filter.incl(1), 0, [], 0);end
    if S.filter.incl(2)>0; EEG = pop_eegfiltnew( EEG, 0, S.filter.incl(2), [], 0);end

    % DOWNSAMPLE
    if S.downsample
        EEG = pop_resample(EEG, S.downsample)
    end

    % EPOCH
    %create epochs if markers exist
    try
        EEG = pop_epoch( EEG, S.epoch.markers, S.epoch.timewin);
    catch
        % if markers don't exist, add markers and epoch
        if S.epoch.addmarker && isempty(EEG.epoch)
            Sr = EEG.srate; % sampling rate of data
            Ndp = Sr*(S.epoch.timewin(2)-S.epoch.timewin(1));% number of data points per epoch
            Tdp = size(EEG.data,2);% total number of data points in file
            Mep = floor(Tdp/Ndp);% max possible number of epochs
            for i = 1:Mep;
                EEG.event(1,i).type = S.epoch.markers{1};
                EEG.event(1,i).latency = Sr*S.epoch.timewin(1)+(i-1)*Ndp+1;
                EEG.event(1,i).urevent = S.epoch.markers{1};
            end
            EEG = pop_epoch( EEG, S.epoch.markers(1), S.epoch.timewin);
        end
    end
    
    sname_ext = 'epoched.set';
    % SEPARATE INTO FILES ACCORDING TO MARKER TYPE AND SAVE
    if ~isempty(S.epoch.separate)
        nfiles = length(S.epoch.separate);
        INEEG =EEG; 
        allmarkers = {INEEG.epoch.eventtype};
        for n = 1:nfiles
            selectmarkers = S.epoch.markers(S.epoch.separate{n});
            markerindex = find(ismember(allmarkers,selectmarkers));
            EEG = pop_select(INEEG,'trial',markerindex);

            % save .set
            sname = [nme '_' S.epoch.separate_suffix{n} '_' sname_ext];
            pop_saveset(EEG,'filename',sname,'filepath',S.setpath); 
        end
    else
        % OR SAVE SINGLE FILE
        sname = [nme '_' sname_ext];
        EEG = pop_saveset(EEG,'filename',sname,'filepath',S.setpath); 
    end
end

% GET FILE LIST
S.filepath = S.setpath;
S.loadext=sname_ext;
S.sessions={};S.blocks={};S.conds={};
S = getfilelist(S);

for f = length(S.filelist)
    file = S.filelist{f};
    EEG = pop_loadset('filename',file,'filepath',S.setpath);

    % LINEAR DETREND
    if S.epoch.detrend
        for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    end

    % REMOVE BASELINE
    if S.epoch.rmbase
        EEG = pop_rmbase( EEG, [S.epoch.timewin(1)*1000    0]);
    end

    EEG = eeg_checkset( EEG );

    % SAVE
    [pth nme ext] = fileparts(file); 
    sname = nme;
    EEG = pop_saveset(EEG,'filename',sname,'filepath',S.setpath); 
end

%% combine data files
if S.combinefiles
    clear OUTEEG
    
    % GET FILE LIST
    S.filepath = S.setpath;
    S.loadext=sname_ext;
    S = getfilelist(S);
    
    for s = 1:length(S.subjects)
        
        % FIND THE FILES FOR THIS SUBJECT
        subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,S.subjects{s})))));
       
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
            EEG = pop_loadset('filename',file,'filepath',S.setpath);
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
        sname_ext = 'combined.set';
        sname = [S.subjects{s} '_' sname_ext];
        EEG = pop_saveset(EEG,'filename',sname,'filepath',S.setpath);
    end
end

% FIND THE NEW DATA FILES
%files = dir(fullfile(S.setpath,['*' S.runsubject '*' sname_ext]));
%files_ana = 1:length(files);

% NOISY TRIAL AND CHANNEL REJECTION USING FIELDTRIP
if ~isempty(S.FTrej) && iscell(S.FTrej)
    
    % GET FILE LIST
    S.filepath = S.setpath;
    S.loadext=sname_ext;
    S.sessions={};S.blocks={};S.conds={};
    S = getfilelist(S);
    
    for f = length(S.filelist)
        file = S.filelist{f};
        EEG = pop_loadset('filename',file,'filepath',S.setpath);

        % strategy - only remove a very small number of very bad trials / chans
        % before ICA - do further cleaning after ICA
        for i = 1:length(S.FTrej)
            EEG = FTrejman(EEG,S.FTrej{i}); 
        end

        % SAVE
        [pth nme ext] = fileparts(file); 
        sname_ext = 'manrej.set';
        sname = [nme '_' sname_ext];
        EEG = pop_saveset(EEG,'filename',sname,'filepath',S.setpath); 
    end
end

% ICA
if S.ICA
    % GET FILE LIST
    S.filepath = S.setpath;
    S.loadext=sname_ext;
    S = getfilelist(S);
    
    for f = length(S.filelist)
        file = S.filelist{f};
        EEG = pop_loadset('filename',file,'filepath',S.setpath);

        %RUN ICA 
        numcomp = numcompeig(EEG);
        EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',numcomp);

        % SAVE
        [pth nme ext] = fileparts(file); 
        sname_ext = 'ICA.set';
        sname = [nme '_' sname_ext];
        EEG = pop_saveset(EEG,'filename',sname,'filepath',S.setpath);
    end
end