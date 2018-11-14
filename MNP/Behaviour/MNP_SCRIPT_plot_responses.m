%% Analysis of response model predictions: fit RT and compare prediction of Choice with BayesOpt predictions
% testing this on VLTv2 cab10 data, response model 10 wins (mostly driven by beta2 parameter)
clear all
close all
dbstop if error

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'VLT'; %
S.version = 2; % version of the NLT/ALT design
% S.path.seq = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\seq']; % unprocessed data in original format
% S.path.raw = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\raw']; % unprocessed data in original format
% S.path.prep = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\processed']; % folder to save processed .set data
% S.path.hgf = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\hgf']; % folder to save processed data
S.path.seq = ['C:\Data\MNP\behaviour\Pilots\seq']; 
S.path.raw = ['C:\Data\MNP\behaviour\Pilots\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\MNP\behaviour\Pilots\processed']; % folder to save processed .set data
S.path.hgf = ['C:\Data\MNP\behaviour\Pilots\hgf']; % folder to save processed data
S.fname.parts = {'prefix','study','subject','block','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.study = {'MNP'};
S.subjdir = 1;
S.select.subjects = {'test'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {['Sequence_' S.expt '_OptionALPL_startblock*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.select.condtype = ''; % specify whether to analyse using AL or PL condnums; blank = use seq.condnum
S.path.datfile = ['C:\Data\MNP\participants\Participant_Data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
run('C:\Data\Matlab\Matlab_files\MNP\MNP_addpaths');

% data import
S.load.prefixes = {'Output','Sequence'};
[S,D]=SCIn_data_import(S);

%PROCESSING OF RESPONSES
S.fitsim=1;
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 1;
S.RT.min = 0.5; % min RT to consider
[S,D_in]=SCIn_data_process(S,D);

% plot responses
[S,D]= MNP_plot_responses(S,D_in)
