%% PREPROCESSING
% 1. filter, epoch, baseline correction, rereference
% 2. reject chans/trials
% 3. perform ICA decomposition for later rejection of eye movement and
% other noise artefact

% Altered from Chris' original 'Preprocess' script 
%%

%clear all
%dbstop if error

% DIRECTORIES CONTAINING DATA
origpath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Set'; % unprocessed data
anapath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Preprocessed'; % folder to save analysed data
design_path = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Raw'; % contains some related information about the study design

% FIND THE DATA FILES
cd(origpath);
files = dir('NTIP_TimEyes*.set');

cd(anapath);

% SET SOME OPTIONS
timebin = [0 0.5]; % for epoching
filterset = [0.5 100]; % FILTER SETTINGS - INCLUSIVE
notch_on = 3;
addpath(genpath('Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Supplementary data'));

files_ana = 1:length(files);
trials_ana = 1; fname_ext = '';

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    % LOAD DATA
    EEG = pop_loadset('filename',orig_file,'filepath',origpath);
    
        %ADD CHANNEL LOCATIONS 
    EEG=pop_chanedit(EEG, 'lookup','C:\\Work\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
    
    % EXCLUDE FC2 BEFORE INTERPOLATION
    chanexcl = [22];   
    EEG= eeg_interp(EEG, chanexcl, 'spherical');
   
    % APPLY NOTCH FILTER: 3 OPTIONS
    if notch_on==1 
        EEG = pop_eegfiltnew(EEG,45,55,[],1); 
        EEG = pop_eegfiltnew(EEG,95,105,[],1); 
    elseif notch_on==2
        EEG.data = rm50Hz(EEG.data,[2 1 3],EEG.srate,20);
    elseif notch_on==3
        %% Parameters for PrepPipeline - performs line denoising and robust re-referencing
        % conventional re-referencing to the common average, prior to
        % cleaning could introduce artefact to otherwise clean channels.
        % Robust averaging as done here gets around this.
        params = struct();
        %params.lineFrequencies = [50, 100, 150, 200];
        params.lineFrequencies = [];
        %params.referenceChannels = 1:64;
        %params.evaluationChannels = 1:64;
        %params.rereferencedChannels = 1:70;
        %params.detrendChannels = 1:70;
        %params.lineNoiseChannels = 1:70;
        params.detrendType = 'high pass';
        params.detrendCutoff = 1;
        params.referenceType = 'robust';
        params.meanEstimateType = 'median';
        params.interpolationOrder = 'post-reference';
        params.keepFiltered = false;
        %params.name = thisName;
        params.ignoreBoundaryEvents = true;
        [EEG, computationTimes] = prepPipeline(EEG, params); 
        %fprintf('Computation times (seconds):\n   %s\n', ...
            %getStructureString(computationTimes));
        EEG = pop_eegfiltnew(EEG,45,55,[],1); 
        EEG = pop_eegfiltnew(EEG,95,105,[],1); 
    end
    
    % FILTER
    if filterset(1)>0; EEG = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
    if filterset(2)>0; EEG = pop_eegfiltnew( EEG, 0, filterset(2), [], 0);end
   
    
    % EPOCH
    %add 2sec epochs by adding markers every second and epochs after that
    
    %add the marker ('M') 
    Sr = EEG.srate; % sampling rate of data
    Ndp = Sr*(timebin(2)-timebin(1));% number of data points per epoch
    Tdp = size(EEG.data,2);% total number of data points in file
    Mep = floor(Tdp/Ndp);% max possible number of epochs
    for i = 1:Mep;
        EEG.event(1,i).type = 'M';
        EEG.event(1,i).latency = (i-1)*Ndp+1;
        EEG.event(1,i).urevent = 'M';
    end
    
    %create 2sec epochs
    EEG = pop_epoch( EEG, {  'M'  }, timebin, 'newname', [C{1} '_' C{2} '_epochs'],'epochinfo', 'yes');
  
    % LINEAR DETREND
    %for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    
    % REMOVE BASELINE
    %EEG = pop_rmbase( EEG, [timebin(1)*100    0]);
    
    % SAVE
    sname = [C{1} '_' C{2} '_epoched.set'];
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
    
end


% NOISY TRIAL AND CHANNEL REJECTION USING FIELDTRIP
for f = files_ana
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    
    % strategy - only remove a very small number of very bad trials / chans
    % before ICA - do further cleaning after ICA
    EEG = FTrejman(EEG,[20 30]); % high freq to identify noise not related to eye movement
    EEG = FTrejman(EEG,[0 5]); % low freq to identify eye movement
    
    %EEG = pop_eegplot(EEG);
    
    sname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
end

% prepare data for ICA
for f = files_ana
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    if notch_on~=3 % only re-reference here only if not already done so using prepPipeline
        %RE-REFERENCE TO COMMON AVERAGE
        EEG = pop_reref( EEG, []);
    end

    %RUN ICA 
    numcomp = numcompeig(EEG);
    EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',numcomp);

    % SAVE
    sname = [C{1} '_' C{2} '_ICA.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
    
end