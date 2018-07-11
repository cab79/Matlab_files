%% Analysis of response model predictions: fit RT and compare prediction of Choice with BayesOpt predictions
% testing this on VLTv2 cab10 data, response model 10 wins (mostly driven by beta2 parameter)
clear all
close all
clear S

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'VLT'; %
S.version = 2; % version of the NLT/ALT design
S.path.seq = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\seq']; % unprocessed data in original format
S.path.raw = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\processed']; % folder to save processed .set data
S.path.hgf = ['C:\Data\MNP\behaviour\Pilots\' S.expt 'v' num2str(S.version) '\hgf']; % folder to save processed data
S.fname.parts = {'prefix','subject','block','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.subjects = {'cab10'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {['Sequence_' S.expt '_OptionAssoc*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\MNP\participants\Participant_Data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
MNP_addpaths

% data import
S.load.prefixes = {'Output','Sequence'};
[S,D]=SCIn_data_import(S);

%PROCESSING OF RESPONSES
S.fitsim=1;
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 0;
S.RT.on = 1;
S.RT.min = 0.5; % min RT to consider
[S,D_in]=SCIn_data_process(S,D);

response_models = [2:12];
for rm = 1:length(response_models)
   
    % Fitting: AL model with RTs
    S.perc_model = 4; 
    S=MNP_perceptual_models(S);
    S.resp_model = response_models(rm); 
    S=MNP_response_models(S);
    S.HGF.plottraj = 0; 
    [D_prep,S] = MNP_prep4HGF(D_in,S,1); % prepare response data for HGF
    % HGF fit to RTs
    S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
    D_fit=HGF_run(D_prep,S,0);
    % HGF Bayes optimal
    S.prc_config = 'GBM_config'; S.obs_config = 'bayes_optimal_binary_config_CAB'; S.nstim=[]; S.bayes_opt=1; 
    D_bo=HGF_run(D_prep,S,0);

    % Simulating: AL model with Choices
    S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
    S.perc_model = 4; 
    S=MNP_perceptual_models(S);
    S.resp_model = 1; % choice model 
    S=MNP_response_models(S);
    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    S.numsimrep = 100; % number of simulations to run per parameter combination
    S.sim=[];
    %S.sim = D_fit.HGFbopars_fit.p_prc;
    [D_fit,S] = HGF_sim(D_fit,S); 
    %S.sim = bopars_bo.p_prc;
    [D_bo,S] = HGF_sim(D_bo,S);

    %compare to ground truth
    buttons = D_prep.Output.Settings.buttonopt;
    target_resp = [2 1]; % which response buttons correspond to simulated response outputs
    actual_choices = [];
    actual_choices(ismember(D_prep.Output.pressbutton,buttons(1)))=1;
    actual_choices(ismember(D_prep.Output.pressbutton,buttons(2)))=2;
    fitted_correct=[];
    bayesopt_correct=[];
    for rep = 1:S.numsimrep
        fitted_choices = [];
        fitted_choices(ismember(D_fit.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(1)))=target_resp(1); % use D_prep.Output.presstrial to ensure correct number of trials
        fitted_choices(ismember(D_fit.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(2)))=target_resp(2);
        bayesopt_choices = [];
        bayesopt_choices(ismember(D_bo.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(1)))=target_resp(1);
        bayesopt_choices(ismember(D_bo.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(2)))=target_resp(2);
        fitted_correct(rep) = sum(fitted_choices == actual_choices)/length(actual_choices);
        bayesopt_correct(rep) = sum(bayesopt_choices == actual_choices)/length(actual_choices);
    end
    means_fitted(rm)=mean(fitted_correct);
    stds_fitted(rm)=std(fitted_correct);
    means_bayesopt(rm)=mean(bayesopt_correct);
    stds_bayesopt(rm)=std(bayesopt_correct);
end

figure
hold on
bar(response_models,means_fitted)
errorbar(response_models,means_fitted,stds_fitted,'.')
plot(xlim,mean(means_bayesopt)*[1 1], 'k--')