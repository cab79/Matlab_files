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
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

%% EEG DATA FILE NAMING
% Name the input files as <study name>_<participant ID>_<sessions name_<block name>_<condition name>
% For example: EEGstudy_P1_S1_B1_C1.set
% Any of the elements can be left out as long as what is there is defined.

%% SET DATA PATHS
clear S
S.path=struct;% clears the field
S.path.main = 'C:\Data\CORE\EEG\ana';
S.path.raw = ''; % unprocessed data in original or BIDS format
S.path.prep = [S.path.main '\prep']; % folder to save processed data
S.path.freq = [S.path.main '\freq']; % folder to save frequency analyses
S.path.tf = [S.path.main '\tf']; % folder to save time-frequency analyses
S.path.erp = [S.path.main '\erp']; % folder to save ERP analyses
S.path.datfile = 'C:\Data\CORE\Participants\Participant_data.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
S.path.locfile = fullfile(eeglab_path,'\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% ERP ANALYSIS 
load(fullfile(S.path.main,'S'))
S.tf=struct;% clears the field
S.path.prep = 'C:\Data\CORE\eeg\ana\prep\cleaned\part4'; % folder to load data from
S.tf.fname.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.tf.fname.ext = {'set'}; % file extension of input data: supports 'vhdr' (Brainvision), 'cnt' (Neuroscan), 'mff' (EGI - requires mffimport2.0 toolbox)
S.tf.study = {};
S.tf.select.groups = {}; % either single/multiple groups, or use * to process all groups in the datfile
S.tf.select.subjects = {}; % either single/multiple subjects, or leave blank to process all subjects in the datfile
S.tf.select.sessions = {}; % Can use * as wildcard
S.tf.select.blocks = {}; % blocks to load (each a separate file). Can use * as wildcard
S.tf.select.conds = {}; % conditions to load (each a separate file). Can use * as wildcard
S.tf.load.suffix = {'4_merged_cleaned'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.tf.select.datatype = 'ERP'; %'Freq','TF','Coh','ERP'
% general settings
S.tf.flipchans = [5:8]; % specific to CORE - fnums to flip
S.tf.epoch.detrend = 0;
S.tf.epoch.basewin = [-0.2 0]; % baseline window
S.tf.epoch.rmbase = 1; % remove baseline prior to frequency/ERP
S.tf.epoch.markers = {'all'}; % marker types to analyse
S.tf.epoch.combinemarkers = 0; %include all events in a single condition.
S.tf.epoch.keeptrials = 'no'; % yes or no
% ERP settings
S.tf.CSD.apply=0; % Apply CSD (leave as 0 unless you know what you are doing!)
S.tf.CSD.montage = 'C:\Data\NTIP\CSDmontage_64.mat';
S=erp_freq_analysis(S)
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% GRAND AVERAGE OF ERP/TF/FREQ DATA
close all
load(fullfile(S.path.main,'S'))
S.ga=struct;% clears the field
S.ga.fname.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.ga.study = {};
S.ga.load.suffix = {'4_merged_cleaned_ERP'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.ga.fname.ext = {'.mat'};% generic file suffix
S.ga.select.groups = {}; % group index, or leave blank to process all 
S.ga.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.ga.select.sessions = {};
S.ga.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.ga.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.ga.select.events = 1;
S.ga.select.chans = {[],[]};
S.ga.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.ga.select.freqs = 0; % select which freq to actally process
S.ga.grand_avg.parts = {}; % which file categories to produce separate grand averages? Blank = one GA for all. Options: {'groups','subjects','sessions','blocks','conds'};
S.ga.grand_avg.weighted = 1; % weighted average according to number of trials
S.ga.grand_avg.outliers = 0; % calculate multivariate outliers (SLOW)
% RUN
S=eeg_grand_average(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable


%% GLOBAL FIELD POWER (GFP) ANALYSIS FOR ERP/TF
close all
load(fullfile(S.path.main,'S.mat'))
S.gfp=struct;% clears the field
S.gfp.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.gfp.select.events = 1;
% Setting for GFP analysis
S.gfp.GAnorm.commonavg = 1; % subtract common average from each electrode
S.gfp.GAnorm.gfp = 1; % divide by the mean GFP over time
S.gfp.GArmbase = 1; % remove grand avg baseline
S.gfp.GFPnorm.type = 'nonorm'; % Options: {'nonorm', 'GFPnorm', 'GFPpriornorm1','GFPpriornorm2','GFPbase','GFPp1norm'}; GFPnorm=1;
S.gfp.epoch.basewin = [-0.2 0]; % baseline window
S.gfp.select.chans = {[],[]};
S.gfp.topXpercent = 100;
S=gfp_analysis(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

