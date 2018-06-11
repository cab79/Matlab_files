clear all
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
restoredefaultpath

%% ADD TOOLBOXES TO MATLAB PATH
eeglab_path = 'C:\Data\Matlab\eeglab14_1_1b'; % path to eeglab toolbox
fieldtrip_path = 'C:\Data\Matlab\fieldtrip-20180320'; % path to fieldtrip toolbox
batchfun_path = 'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_functions';
support_path = 'C:\Data\Matlab\Matlab_files\_generic_eeglab_batch\eeglab_batch_supporting_functions';

addpath(fieldtrip_path); ft_defaults; 
addpath(genpath(eeglab_path)); 
rmpath(genpath(fullfile(eeglab_path,'plugins\Biosig3.3.0'))); % NEW - due to conflict between Fieldtrip and eeglab BIOSIG toolbox
addpath(genpath(batchfun_path));
addpath(genpath(support_path));

%% EEG DATA FILE NAMING
% Name the input files as <study name>_<participant ID>_<sessions name_<block name>_<condition name>
% For example: EEGstudy_P1_S1_B1_C1.set
% Any of the elements can be left out. But all must be separated by underscores.

%% SET DATA PATHS
clear S
S.path=struct;% clears the field
S.path.main = 'C:\Data\temp\Sarah'
S.path.raw = [S.path.main '\Raw']; % unprocessed data in original format
S.path.prep = [S.path.main '\SetC']; % folder to save processed .set data
S.path.freq = [S.path.main '\Frequency']; % folder to save frequency analyses
S.path.tf = [S.path.main '\TF']; % folder to save time-frequency analyses
S.path.erp = [S.path.main '\ERP']; % folder to save ERP analyses
S.path.datfile = [S.path.main '\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
S.path.locfile = fullfile(eeglab_path,'\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% PREPROCESSING
%% SET DATA OPTIONS
load(fullfile(S.path.main,'S'))
S.prep=struct; % clears the field
S.prep.fname.parts = {'study','subject','block','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.prep.study = {'xiipcD2'};
S.prep.load.suffix = {'SetC'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.prep.fname.ext = {'set'};% generic file suffix
S.prep.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.prep.select.sessions = {};
S.prep.select.blocks = {'*'}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.prep.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
%% STEP 3: channel and trial rejection
S.prep.save.suffix = {'manrej'}; % suffix to add to output file name, if needed. 
% S.prep.clean.flatchan.varthresh = 1; % std threshold - less variance than this per trial will be rejected
% S.prep.clean.flatchan.trial_weight = 1;
% S.prep.clean.flatchan.chan_weights = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 20 50 100];
S.prep.clean.FTrej = {[]}; % high freq to identify noise, low freq to identify eye movement
S.prep.startfile = 1; 
S=eeglab_preprocess(S,'rej')
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
%% STEP 4: ICA
S.prep.load.suffix = {'manrej'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.prep.save.suffix = {'ICA'}; % suffix to add to output file name, if needed. 
S.prep.clean.ICA = 1;
S.prep.startfile = 1; 
S=eeglab_preprocess(S,'ICA')
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% ERP ANALYSIS 
load(fullfile(S.path.main,'S'))
S.tf=struct;% clears the field
%S.path.prep = 'C:\EEGLABtraining2018\Data\Preprocessed\separated'; % folder to load data from
S.tf.fname.parts = {'study','subject','session','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.tf.study = {'xiipcD2'};
S.tf.load.suffix = {'SetC'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.tf.fname.ext = {'set'};% generic file suffix
S.tf.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.tf.select.sessions = {'01','02','03'};
S.tf.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.tf.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.tf.select.datatype = 'ERP'; %'Freq','TF','Coh','ERP'
% general settings
S.tf.epoch.basewin = [-3.5 -3]; % baseline window
S.tf.epoch.rmbase = 1; % remove baseline prior to frequency/ERP
S.tf.epoch.markers = {'S  4','S  5','S  8','S  9'}; % marker types to analyse
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
S.ploterp.fname.parts = {'study','subject','block','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.ploterp.study = {'xiipcD2'};
S.ploterp.load.suffix = {'SetC_ERP'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.ploterp.fname.ext = {'mat'};% generic file suffix
S.ploterp.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.tf.select.sessions = {'01','02','03'};
S.ploterp.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.ploterp.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.ploterp.select.events = 1:4;
S.ploterp.select.chans = {[],[]};
S.ploterp.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.ploterp.select.freqs = 0; % select which freq to actally process
S.ploterp.layout = 'acticap-64ch-standard2.mat'; % layout file in Fieldtrip layouts folder
S.ploterp.ylim= [-15,15]; 
S.ploterp.times= {[-3 -2.5],[-0.5 0],[0.2 0.25],[0.3 0.5]}; 
% RUN
S=plot_ERPs(S,'subject');
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% GRAND AVERAGE OF ERP/TF/FREQ DATA
close all
load(fullfile(S.path.main,'S'))
S.ga=struct;% clears the field
S.ga.fname.parts = {'study','subject','session','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'prefix','study','group','subject','session','block','cond','suffix','ext'};
S.ga.study = {'xiipcD2'};
S.ga.load.suffix = {'SetC_ERP'}; % suffix to add to input file name, if needed. Can use * as wildcard
S.ga.fname.ext = {'mat'};% generic file suffix
S.ga.select.groups = {}; % group index, or leave blank to process all 
S.ga.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.ga.select.sessions = {'01','02','03'};
S.ga.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.ga.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.ga.select.events = 1:4;
S.ga.select.chans = {[],[]};
S.ga.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.ga.select.freqs = 0; % select which freq to actally process
S.ga.grand_avg.parts = {}; % which file categories to produce separate grand averages? Blank = one GA for all. Options: {'groups','subjects','sessions','blocks','conds'};
S.ga.grand_avg.weighted = 1; % weighted average according to number of trials
S.ga.grand_avg.outliers = 1; % calculate multivariate outliers (SLOW). 1 = use all data; 2 = calculate for each event separately
% RUN
S=eeg_grand_average(S);
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable
% PLOT
close all
S.ploterp.select.datatype = 'ERP'; % Options: 'TF','ERP','Freq'. NOT YET TESTED ON Freq or TF DATA - WILL PROBABLY FAIL
S.ploterp.select.chans = {[],[]};
S.ploterp.layout = 'acticap-64ch-standard2.mat'; % layout file in Fieldtrip layouts folder
S.ploterp.ylim= [-15,15]; 
S.ploterp.times= {[-2.5 -2],[-0.5 0],[0.2 0.25],[0.3 0.5]}; 
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
