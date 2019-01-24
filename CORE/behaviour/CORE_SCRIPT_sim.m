%% Analysis: Perceptual model fitting

dbstop if error
clear all
%close all
clear S

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'part2'; %
S.path.seq = ['C:\Data\CORE\design']; % unprocessed data in original format
S.path.raw = ['C:\Data\CORE\behaviour\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\CORE\behaviour\processed']; % folder to save processed data
S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; % folder to save processed data
S.path.design = ['C:\Data\CORE\design']; % 
S.fname.parts = {'prefix','subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.groups = {};
S.select.subjects = {'CORE002'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

% unique save name extension
sname = datestr(now,30)
S.sname=sname;

% data import
S.load.prefixes = {'RT','dt'};
S.load.suffix = {'*'};
[S,D]=SCIn_data_import(S);

% response data processing
S.fitsim=2; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 0;
S.RT.on = 0;
S.save.tables = 0;
[S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

% model fitting
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_model=[12]; % 1 3 9 10 11 12
S.resp_model = [15];
S=CORE_perceptual_models(S);
S=CORE_response_models(S);
S.HGF.plottraj = 1; % turn off if doing multiple simulations!

% which parameters?
%S.sim=[]; % specify here? or use generic parameters: S=CORE_sim_parameters(S)
%S=CORE_sim_parameters(S);
S.fitted_hgf = 'CORE_fittedparameters_percmodel128_respmodel15_fractrain0_20181228T193827.mat'; 
%S.fitted_hgf = 'CORE_fittedparameters_percmodel12_bayesopt_20181019T083824.mat';
S=CORE_get_median_params_from_fitted(S,[1:22]); %group 1
D_sim=HGF_run(D_prep,S,1);
S=CORE_get_median_params_from_fitted(S,[23:44]); %group 2
D_sim=HGF_run(D_prep,S,1);
%save(fullfile(S.path.hgf,'sim',['CORE_sim_percmodel' num2str(S.perc_model) '_' S.sname '.mat']), 'D_sim', 'S');