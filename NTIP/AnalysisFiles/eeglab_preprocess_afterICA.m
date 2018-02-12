function S=eeglab_preprocess_afterICA(S)
% PREPROCESS AFTER ICA
% To be run after ICA components have been manually selected for rejection
% via EEGLAB or SASICA, but not yet removed from the data.

% GET FILE LIST
S.filepath = S.setpath;
S = getfilelist(S);

for f = length(S.filelist)
    file = S.filelist{f};
    EEG = pop_loadset('filename',file,'filepath',S.setpath);
    
    % remove selected ICA components (MUST HAVE ALREADY SELECTED THESE MANUALLY)
    if S.ICAremove
        EEG = pop_subcomp( EEG, [], 0); 
    end
    
    % detrend
    if S.detrend
        for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
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
    sname_ext = 'cleaned.set';
    sname = [nme '_' sname_ext];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',S.setpath); 
end

% separate files that were originally combined before ICA
if S.separatefiles
    % GET FILE LIST
    S.filepath = S.setpath;
    S.loadext=['combined_' sname_ext];;
    S = getfilelist(S);

    for f = length(S.filelist)
        file = S.filelist{f};
        INEEG = pop_loadset('filename',file,'filepath',S.setpath);
        [sfiles,ia,ib] = unique({INEEG.epoch.file});
        for s = 1:length(sfiles)
            ind = find(ib==s);
            EEG = pop_select(INEEG,'trial',ind);

            % save .set
            sname = [sfiles{s} '_cleaned.set'];
            pop_saveset(EEG,'filename',sname,'filepath',S.setpath); 
        end
    end
end