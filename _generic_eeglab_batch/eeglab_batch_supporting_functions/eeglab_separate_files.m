
function S=eeglab_separate_files(S)

sname_ext = 'separated';

% separate files that were originally combined before ICA
if S.(S.func).separatefiles.on
    % GET FILE LIST
    S.path.file = fullfile(S.path.prep,'cleaned');
    S = getfilelist(S,'ICA_cleaned');

    loadpath = S.path.file;
    for f = 1:length(S.(S.func).filelist)
        file = S.(S.func).filelist{f};
        INEEG = pop_loadset('filename',file,'filepath',loadpath);
        [sfiles,ia,ib] = unique({INEEG.epoch.file});
        for s = 1:length(sfiles)
            ind = find(ib==s);
            EEGsep1 = pop_select(INEEG,'trial',ind);
            
            sname_ext = 'separated';
            if ~exist(fullfile(S.path.prep,sname_ext),'dir')
                mkdir(fullfile(S.path.prep,sname_ext));
            end
            % SEPARATE INTO FILES ACCORDING TO MARKER TYPE AND SAVE
            if ~isempty(S.(S.func).separatefiles.markerindex)
                nfiles = length(S.(S.func).separatefiles.markerindex);
                INEEGsep1 =EEGsep1; 
                allmarkers = {INEEGsep1.epoch.eventtype};
                for n = 1:nfiles
                    selectmarkers = S.(S.func).epoch.markers(S.(S.func).separatefiles.markindex{n});
                    markerindex = find(ismember(allmarkers,selectmarkers));
                    EEG = pop_select(INEEGsep1,'trial',markerindex);

                    % save .set
                    sname = [sfiles{s} '_' S.(S.func).separatefiles.suffix{n} '.' S.(S.func).fname.ext{:}];
                    pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext)); 
                end
            else
                % save as one file
                sname = [sfiles{s} '_cleaned.' S.(S.func).fname.ext{:}];
                EEG=EEGsep1;
                pop_saveset(EEG,'filename',sname,'filepath',fullfile(S.path.prep,sname_ext));
            end 
        end
         
    end
    
    
end