function S=eeglab_preprocess_afterICA(S)
% PREPROCESS AFTER ICA
% To be run after ICA components have been manually selected for rejection
% via EEGLAB or SASICA, but not yet removed from the data.

% GET FILE LIST
sname_ext = 'ICA';
S.filepath = fullfile(S.setpath,sname_ext);
S = getfilelist(S);

loadpath = fullfile(S.setpath,sname_ext);
for f = 1:length(S.filelist)
    file = S.filelist{f};
    EEG = pop_loadset('filename',file,'filepath',loadpath);
    
    % remove selected ICA components (MUST HAVE ALREADY SELECTED THESE MANUALLY)
    if S.ICAremove
        EEG = pop_subcomp( EEG, [], 0); 
    end
    
    % detrend
    if S.detrend
        for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
    end
    
    % remove baseline
    if S.rmbase
        EEG = pop_rmbase( EEG, [S.basewin(1)*1000 S.basewin(2)*1000]);
    end
    
    % perform final (post-ICA) rejection of bad channels/epochs
    if ~isempty(S.FTrej) && iscell(S.FTrej)
        for i = 1:length(S.FTrej)
            EEG = FTrejman(EEG,S.FTrej{i}); 
        end
    end
    
    % re-reference to the common average
    if S.reref == 1
        % re-reference to the common average
        EEG = pop_reref( EEG, []); 
    elseif S.reref==2
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
   
    % save .set
    [pth nme ext] = fileparts(file); 
    sname_ext = 'cleaned';
    sname = [nme '_' sname_ext '.' S.loadext];
    if ~exist(fullfile(S.setpath,sname_ext),'dir')
        mkdir(fullfile(S.setpath,sname_ext));
    end
    EEG = pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.setpath,sname_ext)); 
    
end

% separate files that were originally combined before ICA
if S.separatefiles
    % GET FILE LIST
    S.filepath = fullfile(S.setpath,sname_ext);
    S = getfilelist(S,sname_ext);

    loadpath = fullfile(S.setpath,sname_ext);
    for f = 1:length(S.filelist)
        file = S.filelist{f};
        INEEG = pop_loadset('filename',file,'filepath',loadpath);
        [sfiles,ia,ib] = unique({INEEG.epoch.file});
        for s = 1:length(sfiles)
            ind = find(ib==s);
            EEGsep1 = pop_select(INEEG,'trial',ind);
            
            sname_ext = 'separated';
            if ~exist(fullfile(S.setpath,sname_ext),'dir')
                mkdir(fullfile(S.setpath,sname_ext));
            end
            % SEPARATE INTO FILES ACCORDING TO MARKER TYPE AND SAVE
            if ~isempty(S.separate)
                nfiles = length(S.separate);
                INEEGsep1 =EEGsep1; 
                allmarkers = {INEEGsep1.epoch.eventtype};
                for n = 1:nfiles
                    selectmarkers = S.epoch.markers(S.separate{n});
                    markerindex = find(ismember(allmarkers,selectmarkers));
                    EEG = pop_select(INEEGsep1,'trial',markerindex);

                    % save .set
                    sname = [sfiles{s} '_' S.separate_suffix{n} '.' S.loadext];
                    pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.setpath,sname_ext)); 
                end
            else
                % save as one file
                sname = [sfiles{s} '_cleaned.' S.loadext];
                EEG=EEGsep1;
                pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.setpath,sname_ext));
            end 
        end
         
    end
    
    
end