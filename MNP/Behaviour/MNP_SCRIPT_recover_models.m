%% Analysis: Perceptual model simulation and recovery
clear all
close all
clear S

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'VLT'; %
S.version = 2; % version of the NLT/ALT design
S.path.seq = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\seq']; % unprocessed data in original format
S.path.raw = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\processed']; % folder to save processed .set data
S.path.hgf = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\hgf']; % folder to save processed .set data
S.fname.parts = {'prefix','subject','block','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.subjects = {'cab15'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {['Sequence_' S.expt '_OptionALPL*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
%S.select.blocks = {['Sequence_' S.expt '_OptionNLT_assoc*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
%S.select.blocks = {['Sequence_' S.expt '_OptionALPL']}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.select.condtype = '';
S.path.datfile = ['C:\Data\MNP\participants\Participant_Data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
MNP_addpaths

% data import
S.load.prefixes = {'Output','Sequence'};
%[S,D]=SCIn_seq_import(S);
[S,D]=SCIn_data_import(S);

% Simulating: AL model with Choices
S.fitsim=2;
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_model = 3; 
S=MNP_perceptual_models(S);
S.resp_model = 1; % choice model 
S=MNP_response_models(S);
S.HGF.plottraj = 0; % turn off if doing multiple simulations!
S.numsimrep = 1; % number of simulations to run per parameter combination
S.sim=[]; % specify here? or use generic parameters: S=MNP_sim_parameters(S)
S=MNP_sim_parameters(S)
[D,S] = MNP_prep4HGF(D,S,2); % prepare response data for HGF
[D,S] = HGF_sim(D,S); 

%PROCESSING OF RESPONSES
S.fitsim=2; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 0;
S.RT.min = 0.5; % min RT to consider
[S,D]=SCIn_data_process(S,D);

% model fitting and recovery
S.perc_model_func = @MNP_perceptual_models;
S.model_sim = S.perc_model; % record previous model simulated
S.perc_model = [3:7]; % select multiple models to test model recovery
[S,D] = HGF_modelrecovery(S,D)
evi=D.evi;
save(fullfile(S.path.prep,['evidence_model' num2str(S.model_sim) '.mat']), 'evi');

% plot responses
%[S,D]= MNP_plot_responses(S,D)

% plot model recovery
if 0
    S.plot_fit_models = 3:7;
    S.plot_sim_models = 3:7;
    HGF_plot_model_recovery(S)
end