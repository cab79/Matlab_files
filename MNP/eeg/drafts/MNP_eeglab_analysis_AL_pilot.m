clear all
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
restoredefaultpath

%% ADD TOOLBOXES TO MATLAB PATH
eeglab_path = 'C:\Data\Matlab\eeglab14_1_1b'; % path to eeglab toolbox
fieldtrip_path = 'C:\Data\Matlab\fieldtrip-20180320'; % path to fieldtrip toolbox
batchfun_path = 'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_functions';
support_path = 'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_supporting_functions';


%% load eeglab (for manual analysis)
if 0
    addpath(eeglab_path); 
    eeglab
end

%% or add paths for script analysis
addpath(fieldtrip_path); ft_defaults; 
addpath(genpath(eeglab_path)); 
rmpath(genpath(fullfile(eeglab_path,'plugins\Biosig3.3.0'))); % due to conflict between Fieldtrip and eeglab BIOSIG toolbox
addpath(genpath(batchfun_path));
addpath(genpath(support_path));

%% EEG DATA FILE NAMING
% Name the input files as <study name>_<participant ID>_<sessions name_<block name>_<condition name>
% For example: EEGstudy_P1_S1_B1_C1.set
% Any of the elements can be left out as long as what is there is defined.

%% SET DATA PATHS
clear S
S.path=struct;% clears the field
S.path.main = 'C:\Data\MNP\eeg\ana';
S.path.raw = 'Q:\MNP\eeg\raw'; % unprocessed data in original or BIDS format
S.path.prep = [S.path.main '\prep']; % folder to save processed data
S.path.freq = [S.path.main '\freq']; % folder to save frequency analyses
S.path.tf = [S.path.main '\tf']; % folder to save time-frequency analyses
S.path.erp = [S.path.main '\erp']; % folder to save ERP analyses
S.path.datfile = 'C:\Data\MNP\participants\participant_data.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
S.path.locfile = fullfile(eeglab_path,'\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% DATA IMPORT
% This identifies files for conversion to EEGLAB format for a
% specific study and sets some parameters for the conversion. It calls the generic function
% "eeglab_import" to do the actual conversion.
load(fullfile(S.path.main,'S'))
S.import=struct;% clears the field
S.import.grpdir = 1; % e.g. 1 if the files are within a group-specific folder with the group's name
S.import.subjdir = 1; % e.g. 1 if the files are within a subject-specific folder with the subject's name
S.import.fname.parts = {'study','subject','block','cond','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.import.fname.partsep = {'_','_','_','',''}; % separator used in filename after each part (must be one per part even if empty)
S.import.fname.ext = {'vhdr'}; % file extension of input data: supports 'vhdr' (Brainvision), 'cnt' (Neuroscan), 'mff' (EGI - requires mffimport2.0 toolbox)
S.import.study = {'MNP'};
S.import.select.groups = {'pilot'}; % either single/multiple groups, or use * to process all groups in the datfile
S.import.select.subjects = {'m001'}; % either single/multiple subjects, or leave blank to process all subjects in the datfile
S.import.select.sessions = {}; % Can use * as wildcard
S.import.select.blocks = {'NLT'}; % blocks to load (each a separate file). Can use * as wildcard
S.import.select.conds = {'ALPL','ALPL_EEG'}; % conditions to load (each a separate file). Can use * as wildcard
S.import.select.chans = {[],[]}; % include (first cell) or exclude (second cell) channels, or leave empty (default all chans)
S.import.chan.addloc = 1; % add channel locations from S.path.locfile; or leave as 0
S.import.load.prefix = {''}; % prefix to add to input file name, if needed. Can use * as wildcard
S.import.load.suffix = {''}; % suffix to add to input file name, if needed. Can use * as wildcard
S.import.save.prefix = {''}; % prefix to add to output file, if needed. 
S.import.save.suffix = {''}; % suffix to add to output file, if needed. 
% RUN
S=eeglab_import(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is re-run

%% PREPROCESSING
%% SET DATA OPTIONS
load(fullfile(S.path.main,'S'))
S.prep=struct; % clears the field
S.prep.loaddir = [S.path.prep '\imported']; % folder to load data from
S.prep.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.prep.fname.partsep = {'_','_','_','_','_','_'}; % separator used in filename after each part (must be one per part even if empty)
S.prep.fname.ext = {'set'}; % file extension of input data: supports 'vhdr' (Brainvision), 'cnt' (Neuroscan), 'mff' (EGI - requires mffimport2.0 toolbox)
S.prep.study = {'MNP'};
S.prep.select.groups = {'pilot'}; % either single/multiple groups, or use * to process all groups in the datfile
S.prep.select.subjects = {'m001'}; % either single/multiple subjects, or leave blank to process all subjects in the datfile
S.prep.select.sessions = {}; % Can use * as wildcard
S.prep.select.blocks = {'NLT'}; % blocks to load (each a separate file). Can use * as wildcard
S.prep.select.conds = {'ALPL','ALPL_EEG'}; % conditions to load (each a separate file). Can use * as wildcard
%% STEP 1: filtering, epoching
S.prep.load.suffix = {''}; % suffix to add to input file name, if needed. Can use * as wildcard
S.prep.save.suffix = {'epoched'}; % suffix to add to output file name, if needed. 
S.prep.cont.timewin = {[],[],[],[]}; % timewindow (s) of data to analyse; one per S.select.conds. Blank = analyse all.
S.prep.cont.downsample = 250;
S.prep.chan.addloc = 0; % add channel locations from this path; or leave as ''
S.prep.chan.switchhead = [];%[1 2]; % changes data channels from headbox 1 to 2 if incorrectly placed during recording
S.prep.chan.interp = [];
S.prep.chan.reref = 0;
S.prep.filter.notch = {[45 55]};
S.prep.filter.incl = [0 95];%[0.1 100]; % FILTER - INCLUSIVE
S.prep.epoch.markers = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15'};
S.prep.epoch.addmarker = 0; % set to 1 to add markers if the above are not present: will use the first marker value
S.prep.epoch.importmarker = 'C:\Data\MNP\behaviour\Pilots\raw\e001\Sequence_MNP_e001_Sequence_NLT_OptionALPL_startblock1_20181003T144514.mat'; % name of file to import markers from; or leave empty
S.prep.epoch.timewin = [-1.5 1]; % peri-marker time window. 
S.prep.epoch.detrend = 1;
S.prep.epoch.rmbase = 1;
S.prep.epoch.basewin = [-1.5 -1]; % baseline window
S.prep.startfile = 1; 
S=eeglab_preprocess(S,'epoch')
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
%% STEP 2: COMBINE FILES
%S.prep.combinefiles.on = 0;
% S.prep.load.suffix = {'epoched'}; % suffix to add to input file name, if needed. Can use * as wildcard
% S.prep.save.suffix = {'combined'}; % suffix to add to output file name, if needed. 
% S.prep.startfile = 1; 
% S=eeglab_preprocess(S,'combine')
% save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
%% STEP 3: channel and trial rejection
S.prep.load.suffix = {'epoched'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.prep.save.suffix = {'manrej'}; % suffix to add to output file name, if needed. 
S.prep.clean.flatchan.varthresh = 1; % std threshold - less variance than this per trial will be rejected
S.prep.clean.flatchan.trial_weight = 1;
S.prep.clean.flatchan.chan_weights = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 20 50 100];
S.prep.clean.FTrej = {[20 90]}; % high freq to identify noise, low freq to identify eye movement
S.prep.startfile = 1; 
S=eeglab_preprocess(S,'rej')
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
%% Separate files into left and right hand
%S.prep.separatefiles.on = 0;
% S.prep.separate.markerindex = {}; % index of markers to separate into files
% S.prep.separatefiles.suffix = {}; % suffixes to new file names: one per cell of S.epoch.separate
% S=eeglab_preprocess(S,'sep')
% save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
%% STEP 4: ICA
S.prep.clean.ICA = 1;
S.prep.load.suffix = {'manrej'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.prep.save.suffix = {'ICA'}; % suffix to add to output file name, if needed. 
S.prep.startfile = 1; 
S=eeglab_preprocess(S,'ICA')
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% PREPROCESSING AFTER ICA
% Only run this part if you have already manually selected ICA components
% for removal and saved that file (with the same name, i.e. ..._ICA.set)
load(fullfile(S.path.main,'S'))
S.path.prep = [S.path.main '\prep']; % folder to save processed data
S.prep2=struct;% clears the field
S.prep2.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.prep2.fname.ext = {'set'};
S.prep2.study = {'MNP'};
S.prep2.load.prefix = {''}; % prefix to add to input file name, if needed. Can use * as wildcard
S.prep2.load.suffix = {'ICA'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.prep2.select.groups = {'pilot'}; % either single/multiple groups, or use * to process all groups in the datfile
S.prep2.select.subjects = {'m001'}; % either a single subject, or leave blank to process all subjects in folder
S.prep2.select.sessions = {};
S.prep2.select.blocks = {'NLT'}; % blocks to load (each a separate file). Can use * as wildcard
S.prep2.select.conds = {'ALPL','ALPL_EEG'}; % conditions to load (each a separate file). Can use * as wildcard
S.prep2.epoch.ICAremove = 1; % remove ICA components (0 if already removed from data, 1 if selected but not removed)
S.prep2.epoch.detrend = 0;
S.prep2.epoch.rmbase = 1;
S.prep2.epoch.basewin = [-1.5 -1]; % baseline window
S.prep2.clean.FTrej.freq = {[]};
S.prep2.clean.FTrej.chan = {[],[]}; % include (first cell) or exclude (second cell) channels, or leave empty (default all chans)
S.prep2.epoch.reref = 1;
S.prep2.separatefiles.on = 0;
S.prep2.separatefiles.markerindex = {}; % index of markers to separate into files
S.prep2.separatefiles.suffix = {}; % suffixes to new files names: one per cell of S.epoch.separate
S.prep2.startfile = 1; 
% RUN
S=eeglab_preprocess_afterICA(S)
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% FREQUENCY / TF ANALYSIS
% load(fullfile(S.path.main,'S'))
% S.tf=struct;% clears the field
% S.path.prep = 'C:\EEGLABtraining2018\Data\Preprocessed\separated'; % folder to load data from
% S.tf.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
% S.tf.study = {'NTIP'};
% S.tf.load.suffix = {'epoched_cleaned'}; % suffix to add to input file name, if needed. Can use * as wildcard
% S.tf.fname.ext = {'set'};% generic file suffix
% S.tf.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
% S.tf.select.sessions = {};
% S.tf.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
% S.tf.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
% S.tf.select.datatype = 'Freq'; %'Freq','TF','Coh','ERP'
% S.tf.select.freq = [6:1:14]; % select frequency range and resolution (if freq analysis). 1D array produces all freqs enteredl 2D array uses ranges, e.g. [0 4; 4 8; 8 13; 13 30; 30 40];
% % general settings
% S.tf.epoch.basewin = [-0.2 0]; % baseline window
% S.tf.epoch.rmbase = 0; % remove baseline prior to frequency/ERP
% S.tf.epoch.markers = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8', 'S  9'}; % marker types to analyse
% S.tf.epoch.combinemarkers = 0; %include all events in a single condition.
% % freq analysis settings
% S.tf.freq.type = 'both'; % induced or both
% S.tf.freq.basenorm = 'relative'; % for TF analysis only: 'absolute', 'relative', 'relchange', 'normchange'
% S.tf.freq.bootrep = 500;% Bootstraps repetitions (for coherence analysis only)
% S.tf.freq.method = 'mtmfft'; % options: 'mtmfft' (Freq),'mtmconvol' (TF), 'wavelet' (TF), 
% S.tf.freq.param=0; % for mtmfft or wavelet, number of cycles. For multitaper, smoothness. If set to zero the time window will be fixed across frequencies to the smallest possible given the freq resolution.
% S.tf.freq.taper = 'hanning'; % 'hanning' or 'dpss'
% S.tf.freq.pad = 'min_pad';
% S.tf.op.mov_avg_trials = 20; % N trials: average EEG power over trials using a moving average window
% S.tf.op.mov_avg_trials_step = 20;
% S=erp_freq_analysis(S)
% save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% ERP ANALYSIS 
load(fullfile(S.path.main,'S'))
S.tf=struct;% clears the field
S.path.prep = 'C:\Data\MNP\eeg\ana\prep\cleaned'; % folder to load data from
S.tf.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.tf.fname.ext = {'set'}; % file extension of input data: supports 'vhdr' (Brainvision), 'cnt' (Neuroscan), 'mff' (EGI - requires mffimport2.0 toolbox)
S.tf.study = {'MNP'};
S.tf.select.groups = {'pilot'}; % either single/multiple groups, or use * to process all groups in the datfile
S.tf.select.subjects = {'m001'}; % either single/multiple subjects, or leave blank to process all subjects in the datfile
S.tf.select.sessions = {}; % Can use * as wildcard
S.tf.select.blocks = {'NLT'}; % blocks to load (each a separate file). Can use * as wildcard
S.tf.select.conds = {'ALPL','ALPL_EEG'}; % conditions to load (each a separate file). Can use * as wildcard
S.tf.load.suffix = {'ICA_cleaned'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.tf.select.datatype = 'ERP'; %'Freq','TF','Coh','ERP'
% general settings
S.tf.epoch.basewin = [-1.5 -1]; % baseline window
S.tf.epoch.rmbase = 1; % remove baseline prior to frequency/ERP
S.tf.epoch.markers = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18','S 19','S 20','S 21','S 22','S 23','S 24'}; % marker types to analyse
S.tf.epoch.combinemarkers = 0; %include all events in a single condition.
% ERP settings
S.tf.CSD.apply=0; % Apply CSD (leave as 0 unless you know what you are doing!)
S.tf.CSD.montage = 'C:\Data\NTIP\CSDmontage_64.mat';
S=erp_freq_analysis(S)
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% PLOT ERPs
close all
load(fullfile(S.path.main,'S'))
S.ploterp=struct;% clears the field
S.ploterp.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.ploterp.study = {'MNP'};
S.ploterp.load.suffix = {'ICA_cleaned_ERP'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.ploterp.fname.ext = {'mat'};% generic file suffix
S.ploterp.select.groups = {'pilot'}; % group index, or leave blank to process all 
S.ploterp.select.subjects = {'m001'}; % either a single subject, or leave blank to process all subjects in folder
S.ploterp.select.sessions = {};
S.ploterp.select.blocks = {'NLT'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.ploterp.select.conds = {'ALPL'}; % conditions to load (each a separate file) - empty means all of them, or not defined
%S.ploterp.select.events = {[1:2:23],[2:2:24]}; % 
%S.ploterp.event_labels = {'low','high'}; % 
%S.ploterp.select.events = {[1 7],[3 5],[2 8],[4 6]}; % 
%S.ploterp.event_labels = {'AL-low-freq','AL-low-rare','AL-high-freq','AL-high-rare'}; % 
%S.ploterp.select.events = {[1 2 7 8],[3 4 5 6]}; % 
%S.ploterp.event_labels = {'AL-freq-cue','AL-rare-cue'}; % 
S.ploterp.select.events = {[9:2:15],[17:2:23],[18:2:24],[10:2:16]}; % 
S.ploterp.event_labels = {'PL-low-freq','PL-low-rare','PL-high-freq','PL-high-rare'}; % 
%S.ploterp.select.events = 1:24; % 
%S.ploterp.event_labels = {}; % 
S.ploterp.select.chans = {[],[]};
S.ploterp.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.ploterp.select.freqs = 0; % select which freq to actally process
S.ploterp.layout = 'acticap-64ch-standard2.mat'; % layout file in Fieldtrip layouts folder
S.ploterp.ylim= [-15,15]; 
S.ploterp.times= {[-1.5 -1],[-1 0],[0.2 0.3],[0.4 0.6]}; 
% RUN
S=plot_ERPs(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% GRAND AVERAGE OF ERP/TF/FREQ DATA
close all
load(fullfile(S.path.main,'S'))
S.ga=struct;% clears the field
S.ga.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.ga.study = {'NTIP'};
S.ga.load.suffix = {'epoched_cleaned_ERP'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.ga.fname.ext = {'.mat'};% generic file suffix
S.ga.select.groups = {'1'}; % group index, or leave blank to process all 
S.ga.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
S.ga.select.sessions = {};
S.ga.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.ga.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.ga.select.events = 1:9;
S.ga.select.chans = {[],[29 30]};
S.ga.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.ga.select.freqs = 10; % select which freq to actally process
S.ga.grand_avg.parts = {}; % which file categories to produce separate grand averages? Blank = one GA for all. Options: {'groups','subjects','sessions','blocks','conds'};
S.ga.grand_avg.weighted = 1; % weighted average according to number of trials
S.ga.grand_avg.outliers = 1; % calculate multivariate outliers (SLOW)
% RUN
S=eeg_grand_average(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
% PLOT
S=plot_ERPs(S,'grandavg');

%% GLOBAL FIELD POWER (GFP) ANALYSIS FOR ERP/TF
close all
load(fullfile(S.path.main,'S.mat'))
S.gfp=struct;% clears the field
S.gfp.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.gfp.select.events = 1:9;
% Setting for GFP analysis
S.gfp.GAnorm.commonavg = 1; % subtract common average from each electrode
S.gfp.GAnorm.gfp = 1; % divide by the mean GFP over time
S.gfp.GArmbase = 1; % remove grand avg baseline
S.gfp.GFPnorm.type = 'nonorm'; % Options: {'nonorm', 'GFPnorm', 'GFPpriornorm1','GFPpriornorm2','GFPbase','GFPp1norm'}; GFPnorm=1;
S.gfp.epoch.basewin = [-0.2 0]; % baseline window
S.gfp.select.chans = {[],[29 30]};
S.gfp.topXpercent = 100;
S=gfp_analysis(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% AVERAGE OVER GROUPS OF CHANNELS: SUBJECTS
close all
load(fullfile(S.path.main,'S.mat'))
S.ca=struct;% clears the field
S.ca.study = {'NTIP'};
S.ca.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.ca.fname.ext = {'mat'};% generic file suffix
S.ca.load.suffix = {'epoched_cleaned_Freq'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.ca.select.datatype = 'Freq'; % Options: 'TF','ERP'. NOT YET TESTED ON TF DATA - WILL PROBABLY FAIL
S.ca.select.groups = {'1'}; % group index, or leave blank to process all 
S.ca.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
S.ca.select.sessions = {};
S.ca.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.ca.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.ca.select.chans = {[],[29 30]}; % ONLY USED IF LOADING A CHANLOCS FILE FOR CHAN SELECTION
S.ca.chanaverage.groups = {[1:20],[21:40],[41:61]}; % either filepath/name to chanlocs file with 'groups' field (one unique value per group); or list of values from 1:length(chanlocs)
% RUN
S=avg_over_chans(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% AVERAGE OVER GROUPS OF CHANNELS: GRAND AVERAGES
% close all
% load(fullfile(S.path.main,'S.mat'))
% S.fname.parts = {'subject'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
% S.fname.ext = 'freq_10.mat';% generic file suffix
% S.select.datatype = 'Freq'; % Options: 'TF','ERP'. NOT YET TESTED ON TF DATA - WILL PROBABLY FAIL
% S.select.groups = {1}; % group index, or leave blank to process all 
% S.select.subjects = {'grandavg'}; % either a single subject, or leave blank to process all subjects in folder
% S.select.sessions = {};
% S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
% S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
% S.select.chans = [29 30]; % ONLY USED IF LOADING A CHANLOCS FILE FOR CHAN SELECTION
% S.chanaverage.groups = {[1:20],[21:40],[41:61]}; % either filepath/name to chanlocs file with 'groups' field (one unique value per group); or list of values from 1:length(chanlocs)
% % RUN
% S=avg_over_chans(S);
% save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% PLOT SINGLE SUBJECT DATA

% 1. PLOT ALL CONDITIONS
close all
load(fullfile(S.path.main,'S'))
S.plot=struct;% clears the field
S.plot.study = {'NTIP'};
S.plot.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.plot.fname.ext = {'mat'};% generic file extension
S.plot.select.datatype = 'Freq'; % Options: 'Freq','TF','ERP'.
% selection of data parts to plot
S.plot.parts = {'conds','chans','freqs'}; % which categories to plot separately? Blank = none (averaged). Options: {'groups','subjects','sessions','blocks','conds','events','freqs','chans'};
S.plot.load.suffix = {'epoched_cleaned_Freq'}; % suffix to add to input file name, if needed. Can use * as wildcard
% select for subplot rows
S.plot.select.groups = {'1'}; % group index, or leave blank to process all 
S.plot.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
S.plot.select.freqs = 10:12; 
% select for subplot columns
S.plot.select.sessions = {};
S.plot.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.plot.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.plot.select.events = 1:9;
S.plot.select.chans = {[],[29 30]};  %first cell: include (default: include all), second cell: exclude (default: exclude none). 
% load and prepare data for plotting
S=eeg_prepare4plot(S);
% plot topographies for each condition
[condlabels,~]=strtok(S.plot.select.conds,'_');
S.plot.type = 'topoplot';
S.plot.title = condlabels;
S.plot.legend= {};
S.plot.scalezero = 0;
S = eeg_plot_func(S);
% all frequencies for each condition
S.plot.type = 'plot_over_freq';
S.plot.title = 'Power at each frequency (no base norm)';
S.plot.legend= condlabels;
S.plot.nele = 5; % num of max power electrodes to avg over
S = eeg_plot_func(S);
% over time for each condition
S.plot.type = 'plot_over_time';
S.plot.title = 'Power over time (no base norm)';
S.plot.legend= condlabels;
S = eeg_plot_func(S);

% 2. PLOT SUBTRACTION BETWEEN CONDITIONS
close all
load(fullfile(S.path.main,'S'))
S.plot=struct;% clears the field
S.plot.study = {'NTIP'};
S.plot.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.plot.fname.ext = {'mat'};% generic file extension
S.plot.select.datatype = 'Freq'; % Options: 'Freq','TF','ERP'.
% selection of data parts to plot
S.plot.parts = {'conds','chans','freqs'}; % which categories to plot separately? Blank = none (averaged). Options: {'groups','subjects','sessions','blocks','conds','events','freqs','chans'};
S.plot.load.suffix = {'epoched_cleaned_Freq'}; % suffix to add to input file name, if needed. Can use * as wildcard
% select for subplot rows
S.plot.select.groups = {'1'}; % group index, or leave blank to process all 
S.plot.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
S.plot.select.freqs = 10:12; 
% select for subplot columns
S.plot.select.sessions = {};
S.plot.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.plot.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.plot.select.events = 1:9;
S.plot.select.chans = {[],[29 30]};  %first cell: include (default: include all), second cell: exclude (default: exclude none). 
% load and prepare data for plotting
S.plot.op.type = 'event';
S.plot.op.operation = 'subtract'; %'subtract', 'relative', 'relchange', 'normchange' or 'db'
S.plot.op.grouping = [1 2; 1 3; 1 4]; % e.g. for subtraction or dividing
S=eeg_prepare4plot(S);
% plot topographies for each condition
S.plot.type = 'topoplot';
S.plot.title = {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
S.plot.legend= {};
S.plot.scalezero = 0;
S = eeg_plot_func(S);
% plot all frequencies for each condition
S.plot.type = 'plot_over_freq';
S.plot.title = 'Subtracted: Entrainment effect at each frequency';
S.plot.legend= {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
S.plot.nele = 5; % num of max power electrodes to avg over
S = eeg_plot_func(S);
%  over time for each condition
S.plot.type = 'plot_over_time';
S.plot.title = 'Subtracted: Entrainment effect over time';
S.plot.legend= {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
S = eeg_plot_func(S);


% 3. PLOT BASELINE NORMALISATION
close all
load(fullfile(S.path.main,'S'))
S.plot=struct;% clears the field
S.plot.study = {'NTIP'};
S.plot.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.plot.fname.ext = {'mat'};% generic file extension
S.plot.select.datatype = 'Freq'; % Options: 'Freq','TF','ERP'.
% selection of data parts to plot
S.plot.parts = {'conds','chans','freqs'}; % which categories to plot separately? Blank = none (averaged). Options: {'groups','subjects','sessions','blocks','conds','events','freqs','chans'};
S.plot.load.suffix = {'epoched_cleaned_Freq'}; % suffix to add to input file name, if needed. Can use * as wildcard
% select for subplot rows
S.plot.select.groups = {'1'}; % group index, or leave blank to process all 
S.plot.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
S.plot.select.freqs = 10:12; 
% select for subplot columns
S.plot.select.sessions = {};
S.plot.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.plot.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.plot.select.events = 1:9;
S.plot.select.chans = {[],[29 30]};  %first cell: include (default: include all), second cell: exclude (default: exclude none). 
% load and prepare data for plotting
S.plot.op.type = 'event';
S.plot.op.operation = 'normchange'; %'subtract', 'relative', 'relchange', 'normchange' or 'db'
S.plot.op.grouping = [1 2; 3 4]; % e.g. for subtraction or dividing
S=eeg_prepare4plot(S);
% Ratio baseline correction: plot topographies for each condition
[condlabels,~]=strtok(S.plot.select.conds,'_');
S.plot.type = 'topoplot';
S.plot.title = condlabels(1:4);
S.plot.legend= {};
S.plot.scalezero = 0;
S = eeg_plot_func(S);
% Ratio baseline correction: plot all frequencies for each condition
S.plot.type = 'plot_over_freq';
S.plot.title = 'Baseline normalisation: Power at each frequency';
S.plot.legend= condlabels(1:4);
S.plot.nele = 5; % num of max power electrodes to avg over
S = eeg_plot_func(S);
% Ratio baseline correction: over time for each condition
S.plot.type = 'plot_over_time';
S.plot.title = 'Baseline normalisation: Power over time';
S.plot.legend= condlabels(1:4);
S = eeg_plot_func(S);

% 4. PLOT BASELINE NORMALISATION AND SUBTRACTION BETWEEN CONDITIONS
close all
load(fullfile(S.path.main,'S'))
S.plot=struct;% clears the field
S.plot.study = {'NTIP'};
S.plot.fname.parts = {'study','subject','block','cond','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.plot.fname.ext = {'mat'};% generic file extension
S.plot.select.datatype = 'Freq'; % Options: 'Freq','TF','ERP'.
% selection of data parts to plot
S.plot.parts = {'conds','chans','freqs'}; % which categories to plot separately? Blank = none (averaged). Options: {'groups','subjects','sessions','blocks','conds','events','freqs','chans'};
S.plot.load.suffix = {'epoched_cleaned_Freq'}; % suffix to add to input file name, if needed. Can use * as wildcard
% select for subplot rows
S.plot.select.groups = {'1'}; % group index, or leave blank to process all 
S.plot.select.subjects = {'P1'}; % either a single subject, or leave blank to process all subjects in folder
S.plot.select.freqs = 10:12; 
% select for subplot columns
S.plot.select.sessions = {};
S.plot.select.blocks = {'ECA'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.plot.select.conds = {'1','1O','10','10O'}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.plot.select.events = 1:9;
S.plot.select.chans = {[],[29 30]};  %first cell: include (default: include all), second cell: exclude (default: exclude none). 
% load and prepare data for plotting
S.plot.op(1).type = 'event';
S.plot.op(1).operation = 'normchange'; %'subtract', 'relative', 'relchange', 'normchange' or 'db'
S.plot.op(1).grouping = [1 2; 3 4]; % e.g. for subtraction or dividing
S.plot.op(2).type = 'event';
S.plot.op(2).operation = 'subtract'; %'subtract', 'relative', 'relchange', 'normchange' or 'db'
S.plot.op(2).grouping = [1 2]; % e.g. for subtraction or dividing
S=eeg_prepare4plot(S);
% Ratio baseline correction: plot topographies for each condition
[condlabels,~]=strtok(S.plot.select.conds,'_');
S.plot.type = 'topoplot';
S.plot.title = {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
S.plot.legend= {};
S.plot.scalezero = 0;
S = eeg_plot_func(S);
% Ratio baseline correction: plot all frequencies for each condition
S.plot.type = 'plot_over_freq';
S.plot.title = 'Subtracted: Entrainment effect at each frequency';
S.plot.legend= {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'}
S.plot.nele = 5; % num of max power electrodes to avg over
S = eeg_plot_func(S);
% Ratio baseline correction: over time for each condition
S.plot.type = 'plot_over_time';
S.plot.title = 'Subtracted: Entrainment effect over time';
S.plot.legend= {'entrain 10Hz-1Hz','entrain 10Hz-1Hz & oddball','mono'};
S = eeg_plot_func(S);

%% Grand average plots
S.gaplot ='freqtopo'; % plot grandavg. Options: 'plottopo' (ERP/TF), 'timtopo' (ERP/TF), 'freqtopo' (Freq)
% TBC

%% EXPORT DATA TO EXCEL FOR FURTHER PLOTTING AND STATISTICS
