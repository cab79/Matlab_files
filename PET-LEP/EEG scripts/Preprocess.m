%% PREPROCESSING
% 1. filtering, epoching, baseline correction, rereference
% 2. rejects chans/trials only with high frequency noise
% 3. reject ICA components related to eye movement (and 50Hz if no notch applied) but not that related to individual channel noise 
% 4. reject remaining chans/trials outside normal range
%%


clear all
filepath = 'C:\Data\PET-LEP\Preprocessed';
cd(filepath);
files = dir('*orig.set');
load chanlocs
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.

timebin= [-3.5 2]; % for epoching, TSOT(4)
basebin = [-3.5 -3];
stimtypes = {'S  1','S  2'};
%ISIs = [1, 0.4];
filterset = [0 50];
notch_on = 1;
addpath(genpath('M:\Matlab\Matlab_files\CORE\Supporting functions'));

for f = 1:length(files)
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    EEG = pop_loadset('filename',orig_file,'filepath',filepath);
    EEG.chanlocs=chanlocs;
    
    % FILTER
    if notch_on==1 
        EEG = pop_eegfiltnew(EEG,45,55,[],1); 
        EEG = pop_eegfiltnew(EEG,95,105,[],1); 
    elseif notch_on==2
        EEG.data = rm50Hz(EEG.data,[2 1 3],EEG.srate,20);
        %pause
        close all
    end;
    if filterset(1)>0; EEG = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
    if filterset(2)>0; EEG = pop_eegfiltnew( EEG, 0, filterset(2), [], 0);end
    
    % EPOCH
    EEG = pop_epoch( EEG, stimtypes, timebin, 'newname', [C{1} '_' C{2} '_epochs'], 'epochinfo', 'yes');
        
    % RE-REFERENCE
    EEG = pop_reref( EEG, []);
    
    % LINEAR DETREND
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    
    % REMOVE BASELINE
    EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
    
    sname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
    
end

for f = 1:length(files)
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_loadset('filename',lname,'filepath',filepath);
    
    EEG = FTrejman(EEG,[20 40]); % high freq to identify noise not related to eye movement
    
    sname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
end

for f = 1:length(files)
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_loadset('filename',lname,'filepath',filepath);
    
    %REJECT CHANNELS
   % [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');
    
    %REJECT EPOCHS
   % EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',12);
    
    %REJECT CHANNELS
   % [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');
    
   EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',30);
   sname = [C{1} '_' C{2} '_ICA.set'];
   EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

   clear EEG
end
eeglab
