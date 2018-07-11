%% Analysis: Perceptual model simulation and recovery
% only need to run on one subject/design, as all subjects would provide
% same results (all based on simulation of a single common design)

dbstop if error
clear all
close all
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
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
CORE_addpaths

% data import
S.load.prefixes = {'RT','dt'};
S.load.suffix = {'*'};
[S,D]=SCIn_data_import(S);

% response data processing
S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 1;
S.RT.min = 0.2; % min RT to consider
S.save.tables = 0;
[S,D]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

%% Predicting actual responses from models: Produce N (e.g. 100) simulations for all subjects
% PL model with Choices
S.fitsim=2;
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
perc_models=1:6
sname = datestr(now,30);
sim_correct = [];
for pm=1:length(perc_models)
    S.perc_model = perc_models(pm); 
    S=CORE_perceptual_models(S);
    S.resp_model = 1; % choice model 
    S=CORE_response_models(S);
    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    S.numsimrep = 100; % number of simulations to run per parameter combination
    S.sim=[]; % specify here? or use generic parameters: S=CORE_sim_parameters(S)
    S=CORE_sim_parameters(S)
    
    for d = 1:length(D)
        
        % simulate
        [D_sim(d),S] = HGF_sim(D(d),S); 
        
        % compare to actual responses
        actual_choices = D(d).HGF.y(:,1)
        for rep = 1:S.numsimrep
            sim_choices = D_sim(d).HGF(rep).sim.y(:,1);
            sim_correct(d,rep) = sum(sim_choices == actual_choices)/length(actual_choices);
        end
    end
    means_sim(pm)=mean(mean(sim_correct,2),1);
    stds_sim(pm)=std(mean(sim_correct,2));
end

figure
hold on
bar(perc_models,means_sim)
errorbar(perc_models,means_sim,stds_sim,'.')

%% Model recovery: Produce N (e.g. 10) simulations for one subject only and fit these
% PL model with Choices
S.fitsim=2;
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
perc_models=1:6;
sname = datestr(now,30);
D=D(1); % select one subject (assumes all subjects have same design)
for pm=1:length(perc_models)
    S.perc_model = perc_models(pm); 
    S=CORE_perceptual_models(S);
    S.resp_model = 1; % choice model 
    S=CORE_response_models(S);
    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    S.numsimrep = 10; % number of simulations to run per parameter combination
    S.sim=[]; % specify here? or use generic parameters: S=CORE_sim_parameters(S)
    S=CORE_sim_parameters(S);
    [D_sim,S] = HGF_sim(D,S); 

    %PROCESSING OF RESPONSES
    S.fitsim=2; % is the analysis of recorded data (1) or simulated data (2)?
    S.meansim=1; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
    S.accuracy.on = 1;
    S.RT.on = 0;
    S.save.tables = 0;
    S.RT.min = 0.2; % min RT to consider
    [S,D_sim]=CORE_data_process(S,D_sim);

    % model fitting and recovery
    S.perc_model_func = @CORE_perceptual_models;
    S.model_sim = S.perc_model; % record previous model simulated
    S.perc_model = [2:6]; % select multiple models to test model recovery
    S.HGF.plottraj = 0; 
    [S,D_sim,evi] = HGF_modelrecovery(S,D_sim)
    save(fullfile(S.path.prep,['evidence_model' num2str(S.model_sim) '_' sname '.mat']), 'evi');
    
end

% plot model recovery
S.plot_sim_models = 1:6;
S.plot_fit_models = 2:6;
HGF_plot_model_recovery(S)