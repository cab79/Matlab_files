% PREPROCESS PART 2
% To be run after ICA components have been manually selected for rejection
% via EEGLAB or SASICA, but not yet removed from the data.

clear all
% path to ICA processed data
filepath = 'C:\Data\NTIP\Preprocessed';
cd(filepath);
% file suffix common to all files needing furter processing
files = dir('*_combined.set');
% path to chanlocs file
%load('Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\chanlocs.mat');

% baseline range
basebin = [-0.2 0];
% stimulus markers to include in analysis
stimtypes = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8' 'S  9'};

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
    % EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    
    % detrend, remove baseline and re-reference to the common average
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
    EEG = pop_reref( EEG, []);
    
    % perform final (poast-ICA) rejection of bad channels/epochs
    EEG = FTrejman(EEG,[0 0]);
    
    % ALTERNATIVE DETREND - don't use
    %for i = 1:EEG.trials
    %    EEG.data(:,:,i) = detrend4(EEG.data(:,:,i)',1,250,250)';
    %end
   
    % save .set
    sname = [C{1} '_cleaned.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

end

% separate files
for f = files_ana
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_cleaned.set'];
    INEEG = pop_loadset('filename',lname,'filepath',filepath);
    [sfiles,ia,ib] = unique({INEEG.epoch.file});
    for s = 1:length(sfiles)
        ind = find(ib==s);
        EEG = pop_select(INEEG,'trial',ind);
        
        % save .set
        sname = [sfiles{s} '_cleaned.set'];
        pop_saveset(EEG,'filename',sname,'filepath',filepath); 
    end
end