%% Analysis: Perceptual model simulation and recovery
clear all
close all
clear S

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'VLT'; %
S.version = 2; % version of the NLT/ALT design
S.path.seq = ['C:\Data\MNP\Pilots\' S.expt 'v' num2str(S.version) '\seq']; % unprocessed data in original format
S.path.raw = ['C:\Data\MNP\Pilots\' S.expt 'v' num2str(S.version) '\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\MNP\Pilots\' S.expt 'v' num2str(S.version) '\processed']; % folder to save processed .set data
S.fname.parts = {'prefix','subject','block','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.subjects = {'cab10'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {['Sequence_' S.expt '_OptionAssoc*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
%S.select.blocks = {['Sequence_' S.expt '_OptionNLT_assoc*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
%S.select.blocks = {['Sequence_' S.expt '_OptionALPL']}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\MNP\Pilots\Participant_Data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
MNP_addpaths

% data import
S.load.prefixes = {'Output','Sequence'};
%[S,D]=SCIn_seq_import(S);
[S,D]=SCIn_data_import(S);

% Simulating: AL model with Choices
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_model = 4; 
S=MNP_perceptual_models(S);
S.resp_model = 1; % choice model 
S=MNP_response_models(S);
S.HGF.plottraj = 0; % turn off if doing multiple simulations!
S.numsimrep = 10; % number of simulations to run per parameter combination
S.sim_param = []; % use default
[D,S,sim] = MNP_sim_HGF(D,S); 

%PROCESSING OF RESPONSES
S.fitsim=2; % is the analysis of recorded data (1) or simulated data (2)?
S.meanD=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 0;
S.RT.min = 0.5; % min RT to consider
[S,D]=SCIn_data_process(S,D);

% model fitting and recovery
S.model_sim = S.perc_model; % record previous model simulated
S.perc_model = [3:7]; % select multiple models to test model recovery
[S,D,fit,evi] = MNP_HGF_modelrecovery(S,D)

% plot responses
%[S,D]= MNP_plot_responses(S,D)

% plot model recovery
S.plot_sim_models = 3:7;
MNP_plot_model_recovery(S)

%% Analysis: fit RT and compare prediction of Choice with BayesOpt predictions
% testing this on VLTv2 cab10 data, response model 10 wins (mostly driven by beta2 parameter)
clear all
close all
clear S

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'VLT'; %
S.version = 2; % version of the NLT/ALT design
S.path.seq = ['C:\Data\MNP\Pilots\' S.expt 'v' num2str(S.version) '\seq']; % unprocessed data in original format
S.path.raw = ['C:\Data\MNP\Pilots\' S.expt 'v' num2str(S.version) '\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\MNP\Pilots\' S.expt 'v' num2str(S.version) '\processed']; % folder to save processed .set data
S.fname.parts = {'prefix','subject','block','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.subjects = {'cab10'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {['Sequence_' S.expt '_OptionAssoc*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\MNP\Pilots\Participant_Data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
MNP_addpaths

% data import
S.load.prefixes = {'Output','Sequence'};
[S,D]=SCIn_data_import(S);

%PROCESSING OF RESPONSES
S.fitsim=1;
S.meanD=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 0;
S.RT.on = 1;
S.RT.min = 0.5; % min RT to consider
[S,D]=SCIn_data_process(S,D);

response_models = [2:12];
for rm = 1:length(response_models)
    
    D(rm)=D(1);
   
    % Fitting: AL model with RTs
    S.perc_model = 4; 
    S=MNP_perceptual_models(S);
    S.resp_model = response_models(rm); 
    S=MNP_response_models(S);
    S.HGF.plottraj = 0; 
    % HGF fit to RTs
    S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
    bopars_fit=MNP_HGF(D(rm),S,[]);
    % HGF Bayes optimal
    S.prc_config = 'GBM_config'; S.obs_config = 'bayes_optimal_binary_config_CAB'; S.nstim=[]; S.bayes_opt=1; 
    bopars_bo=MNP_HGF(D(rm),S,[]);

    % Simulating: AL model with Choices
    S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
    S.perc_model = 4; 
    S=MNP_perceptual_models(S);
    S.resp_model = 1; % choice model 
    S=MNP_response_models(S);
    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    S.numsimrep = 100; % number of simulations to run per parameter combination
    S.sim_param = bopars_fit.p_prc;
    [Dsim_fitted,S,sim] = MNP_sim_HGF(D(rm),S); 
    S.sim_param = bopars_bo.p_prc;
    [Dsim_bayesopt,S,sim] = MNP_sim_HGF(D(rm),S);

    %compare to ground truth
    buttons = D(1).Output.Settings.buttonopt;
    target_resp = [2 1]; % which response buttons correspond to simulated response outputs
    actual_choices = [];
    actual_choices(ismember(D(rm).Output.pressbutton,buttons(1)))=1;
    actual_choices(ismember(D(rm).Output.pressbutton,buttons(2)))=2;
    for rep = 1:S.numsimrep
        fitted_choices = [];
        fitted_choices(ismember(Dsim_fitted(rep).Output.pressbutton(D(rm).Output.presstrial),buttons(1)))=target_resp(1);
        fitted_choices(ismember(Dsim_fitted(rep).Output.pressbutton(D(rm).Output.presstrial),buttons(2)))=target_resp(2);
        bayesopt_choices = [];
        bayesopt_choices(ismember(Dsim_bayesopt(rep).Output.pressbutton(D(rm).Output.presstrial),buttons(1)))=target_resp(1);
        bayesopt_choices(ismember(Dsim_bayesopt(rep).Output.pressbutton(D(rm).Output.presstrial),buttons(2)))=target_resp(2);
        fitted_correct = sum(fitted_choices == actual_choices)/length(actual_choices);
        bayesopt_correct = sum(bayesopt_choices == actual_choices)/length(actual_choices);
        D(rm).Processed.fitted_correct.data(rep) = fitted_correct;
        D(rm).Processed.bayesopt_correct.data(rep) = bayesopt_correct;
    end
    D(rm).Processed.fitted_correct.mean = mean(D(rm).Processed.fitted_correct.data);
    D(rm).Processed.bayesopt_correct.mean = mean(D(rm).Processed.bayesopt_correct.data);
    D(rm).Processed.fitted_correct.std = std(D(rm).Processed.fitted_correct.data);
    D(rm).Processed.bayesopt_correct.std = std(D(rm).Processed.bayesopt_correct.data);
end
for rm = 1:length(response_models)
    means_fitted(rm)=D(rm).Processed.fitted_correct.mean;
    stds_fitted(rm)=D(rm).Processed.fitted_correct.std;
    means_bayesopt(rm)=D(rm).Processed.bayesopt_correct.mean;
    stds_bayesopt(rm)=D(rm).Processed.bayesopt_correct.std;
end

figure
hold on
bar(response_models,means_fitted)
errorbar(response_models,means_fitted,stds_fitted,'.')
plot(xlim,mean(means_bayesopt)*[1 1], 'k--')