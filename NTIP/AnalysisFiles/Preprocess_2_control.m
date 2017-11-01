% PREPROCESS PART 2 CONTROL
% To be run after ICA components have been manually selected for rejection
% via EEGLAB or SASICA, but not yet removed from the data.

clear all
% path to ICA processed data
filepath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Preprocessed';
cd(filepath);
% file suffix common to all files needing furter processing
files = dir('NTIP_TimEyesClosed_ICA.set');
% path to chanlocs file
%load('Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\chanlocs.mat');

ALLEEG_save = 1; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

% baseline range
basebin = [-0.2 -0.1];
% stimulus markers to include in analysis
stimtypes = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8', 'S  9'};

% range of files within 'files' to include in this run
files_ana = 1:length(files);

% run loop
for f = files_ana
    % extract parts of filename for use later, e.g. saving
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    % remove selected ICA components (MUST HAVE ALREADY SELECTED THESE MANUALLY)
    % (take out this line if ICA component already manually removed!)
    EEG = pop_subcomp( EEG, [], 0); 
    
    % interpolate any missing channels using chanlocs
    %EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    
    % detrend, remove baseline and re-reference to the common average
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    %EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]); %baseline
    %correction on resting state data wouldn't work
    EEG = pop_reref( EEG, []);
    
    % perform final (poast-ICA) rejection of bad channels/epochs
    EEG = FTrejman(EEG,[0 0]);
    
    % ALTERNATIVE DETREND - don't use
    %for i = 1:EEG.trials
    %    EEG.data(:,:,i) = detrend4(EEG.data(:,:,i)',1,250,250)';
    %end
   
    % save .set
    sname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

    % Also save each condition ('stimtype') in a separate EEG structure (EEGLAB
    % format), all inside one big structure called ALLEEG
    % This can be used for checking data quality later (use "Plot_ALLEEG.m", but is not used for
    % any further processing.
    if ALLEEG_save
        EEGall=EEG;
        if ALLEEG_save==1
            ALLEEG=struct;
        end
        for st = 1:length(stimtypes)
            if any(strcmp({EEGall.event.type},stimtypes{st}))
                EEG = pop_selectevent(EEGall,'type',stimtypes{st});
                if ALLEEG_save==1
                    if st==1
                        ALLEEG = EEG; 
                    else 
                        ALLEEG(st)=EEG;
                    end
                end
            end
        end
        if ALLEEG_save==1
            sname = [C{1} '_' C{2} '_ALLEEG.mat'];
            save(sname,'ALLEEG');
        end

        if ALLEEG_save==2
            if f==files_ana(1)
                ALLEEG = EEG; 
            else 
                ALLEEG(length(ALLEEG)+1)=EEG;
            end
        end
    end
end
if ALLEEG_save==2
    sname = ['Allfiles_nochange_ALLEEG.mat'];
    save(sname,'ALLEEG');
end