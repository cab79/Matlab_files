clear all
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
restoredefaultpath

%% 1. ADD TOOLBOXES TO MATLAB PATH
eeglab_path = 'C:\Data\Matlab\eeglab14_1_1b'; % path to eeglab toolbox
fieldtrip_path = 'C:\Data\Matlab\fieldtrip-20180320'; % path to fieldtrip toolbox
batchfun_path = 'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch';
support_path = 'C:\Data\Matlab\eeglab_batch_supporting_functions';

addpath(fieldtrip_path); ft_defaults; 
addpath(genpath(eeglab_path)); 
rmpath(genpath(fullfile(eeglab_path,'plugins\Biosig3.3.0'))); % NEW - due to conflict between Fieldtrip and eeglab BIOSIG toolbox
addpath(genpath(batchfun_path));
addpath(genpath(support_path));
%% 2. FOLDER AND FILENAME DEFINITIONS

% FILE NAMING
% Name the input files as <study name>_<participant ID>_<sessions name_<block name>_<condition name>
% For example: EEGstudy_P1_S1_B1_C1.set
% Any of the elements can be left out. But all must be separated by underscores.

clear S
S.rawpath = ''; % unprocessed data in original format
S.setpath = 'C:\Data\PET-LEP\Preprocessed'; % folder to save processed .set data
S.freqpath = 'C:\Data\PET-LEP\Frequency'; % folder to save processed .set data
S.erppath = 'C:\Data\PET-LEP\ERP'; % folder to save processed .set data
S.fnameparts = {'subject'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.subjects = {'P35'}; % either a single subject, or leave blank to process all subjects in folder
S.sessions = {};
S.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.datfile = 'C:\Data\PET-LEP\Participant_data.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.setpath,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% 3. DATA IMPORT
% This identifies files for conversion to EEGLAB format for a
% specific study and sets some parameters for the conversion. It calls the generic function
% "import_eeglab" to do the actual conversion.

% SETTINGS
if 0
load(fullfile(S.setpath,'S'))
S.loadext = 'vhdr'; % file extension of input data: supports 'vhdr' (Brainvision), 'cnt' (Neuroscan), 'mff' (EGI - requires mffimport2.0 toolbox)
S.saveprefix = ''; % prefix to add to output file, if needed
S.savesuffix = ''; % suffix to add to output file, if needed
S.chan.excl = []; % exclude channels, or leave empty as []
S.chan.addloc = fullfile(eeglab_path,'plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp'); % add channel locations from this path; or leave as ''
% RUN
S=eeglab_import(S);
save(fullfile(S.setpath,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
end

%% 2. PREPROCESSING
% SET PREPROCESSING OPTIONS
load(fullfile(S.setpath,'S'))
S.setpath = 'C:\Data\PET-LEP\Preprocessed'; % folder to save processed .set data
S.loadext = 'set';
S.cont.timewin = {}; % timewindow (s) of data to analyse; one per S.conds. Blank = analyse all.
S.downsample = 125;
S.chan.addloc = fullfile(eeglab_path,'plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp'); % add channel locations from this path; or leave as ''
S.chan.interp = [];
S.chan.reref = 1;
S.filter.notch = {};
S.filter.incl = [0 45];%[0.1 100]; % FILTER - INCLUSIVE
S.epoch.markers = {'S  1' 'S  2'};
S.epoch.addmarker = 0; % set to 1 to add markers if the above are not present: will use the first marker value
S.epoch.timewin = [-4 2]; % peri-marker time window
S.epoch.detrend = 1;
S.epoch.rmbase = 1;
S.epoch.separate = {[]}; % index of markers to separate into files
S.epoch.separate_suffix = {}; % suffixes to new file names: one per cell of S.epoch.separate
S.FTrej = {[0 2],[20 40]}; % high freq to identify noise, low freq to identify eye movement
S.ICA = 1;
S.combinefiles = 0;
% RUN
S=eeglab_preprocess(S)
save(fullfile(S.setpath,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% 3. PREPROCESSING AFTER ICA
% Only run this part if you have already manually selected ICA components
% for removal and saved that file (with the same name, i.e. ..._ICA.set)
load(fullfile(S.setpath,'S'))
S.loadext = 'set';
S.subjects = {'S18','P06','P08','P35'}; % either a single subject, or leave blank to process all subjects in folder
%S.subjects = {};
S.ICAremove = 1; % remove ICA components (0 if already removed from data, 1 if selected but not removed)
S.detrend = 1;
S.rmbase = 1;
S.basewin = [-4 -3.5]; % baseline window
S.FTrej = {[]};
S.reref = 1;
S.separatefiles = 0;
S.separate = {}; % index of markers to separate into files
S.separate_suffix = {}; % suffixes to new files names: one per cell of S.epoch.separate
% RUN
S=eeglab_preprocess_afterICA(S)
save(fullfile(S.setpath,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% 4. FREQUENCY ANALYSIS
% directories (if changed)
S.setpath = 'C:\Data\PET-LEP\Preprocessed\cleaned'; % folder to save processed .set data
% data
S.loadext = 'set';% generic file suffix
S.subjects = {'S18','P06','P08','P35'}; % either a single subject, or leave blank to process all subjects in folder
% general settings
S.basewin = [-4 -3.5]; % baseline window
S.rmbase = 0; % remove baseline prior to frequency/ERP
S.markers = {'S  1','S  2'}; % marker types to analyse
S.combinemarkers = 0; %include all events in a single condition.
% analysis type
S.analysistype = 'TF'; %'Freq','TF','Coh','ERP'
% freq analysis settings
S.freqtype = 'both'; % induced or both
S.baselinetype = ''; % for TF analysis only: 'absolute', 'relative', 'relchange', 'normchange'
S.freqsrange = [3:40]; % select frequency range and resolution (if freq analysis). 1D array produces all freqs enteredl 2D array uses ranges, e.g. [0 4; 4 8; 8 13; 13 30; 30 40];
S.mov_avg_trials = 0; % N trials: average EEG power over trials using a moving average window
S.mov_avg_trials_step = 0;
S.bootrep = 0;% Bootstraps repetitions (for coherence analysis only)
S.freq_method = 'mtmconvol'; % options: 'mtmfft' (Freq),'mtmconvol' (TF), 'wavelet' (TF), 
S.freq_param=5; % for mtmfft or wavelet, number of cycles. For multitaper, smoothness. If set to zero the time window will be fixed across frequencies to the smallest possible given the freq resolution.
S.freq_taper = 'hanning'; % 'hanning' or 'dpss'
S.pad = 'min_pad'; %'min_pad', 'nextpow2', 'maxperlen' or a number in seconds
% ERP settings
S.CSD.apply=0; % Apply CSD (leave as 0 unless you know what you are doing!)
S.CSD.montage = 'C:\Data\NTIP\CSDmontage_64.mat';
S=erp_freq_analysis(S)
%% 4. ERP ANALYSIS
% analysis type
S.analysistype = 'ERP'; %'Freq','TF','Coh','ERP'
S=erp_freq_analysis(S)
save(fullfile(S.setpath,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
