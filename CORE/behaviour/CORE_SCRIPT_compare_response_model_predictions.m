%% Analysis of response model predictions: fit RT and compare prediction of Choice with BayesOpt predictions
% testing this on VLTv2 cab10 data, response model 10 wins (mostly driven by beta2 parameter)
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
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
CORE_addpaths

% data import
S.load.prefixes = {'RT','dt'};
S.load.suffix = {'*'};
[S,D]=SCIn_data_import(S);

% response data processing
S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 0;
S.RT.on = 1;
S.RT.min = 0.2; % min RT to consider
S.save.tables = 0;
[S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

% responses
% buttons = D_prep.Output.Settings.buttonopt;
% target_resp = [2 1]; % which response buttons correspond to simulated response outputs
%actual_choices = D_prep.HGF.y(:,1);
% actual_choices(ismember(D_prep.Output.pressbutton,buttons(1)))=1;
% actual_choices(ismember(D_prep.Output.pressbutton,buttons(2)))=2;

sname = datestr(now,30);
response_models = [2:12];
for rm = 4:length(response_models)
   
    % Fitting: model with RTs
    S.perc_model = 6; 
    S=CORE_perceptual_models(S);
    S.resp_model = response_models(rm); 
    S=CORE_response_models(S);
    S.HGF.plottraj = 0; 
    % HGF fit to RTs
    S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
    D_fit=HGF_run(D_prep,S,0);
    % % HGF Bayes optimal
    % S.prc_config = 'GBM_config'; S.obs_config = 'bayes_optimal_binary_config_CAB'; S.nstim=[]; S.bayes_opt=1; 
    % D_bo=HGF_run(D_prep,S,0);

    % Simulating: AL model with Choices
    S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
    S=CORE_perceptual_models(S);
    S.resp_model = 1; % choice model 
    S=CORE_response_models(S);
    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    S.numsimrep = 100; % number of simulations to run per parameter combination
    S.sim=[];
    [D_fit,S] = HGF_sim(D_fit,S); 
    %[D_bo,S] = HGF_sim(D_bo,S);
    
    for d = 1:length(D)
        %compare to ground truth
        %fitted_correct=[];
    %     bayesopt_correct=[];
        actual_choices = D_prep(d).HGF.y(:,1);
        for rep = 1:S.numsimrep
            fitted_choices = D_fit(d).HGF(rep).sim.y(:,1);
    %         fitted_choices(ismember(D_fit.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(1)))=target_resp(1); % use D_prep.Output.presstrial to ensure correct number of trials
    %         fitted_choices(ismember(D_fit.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(2)))=target_resp(2);
            fitted_correct(rep) = sum(fitted_choices == actual_choices)/length(actual_choices);
    %         bayesopt_choices = [];
    %         bayesopt_choices(ismember(D_bo.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(1)))=target_resp(1);
    %         bayesopt_choices(ismember(D_bo.Output(rep).pressbutton(D_prep.Output.presstrial),buttons(2)))=target_resp(2);
    %         bayesopt_correct(rep) = sum(bayesopt_choices == actual_choices)/length(actual_choices);
        end
        means_fitted(d,rm)=mean(fitted_correct);
        stds_fitted(d,rm)=std(fitted_correct);
    %     means_bayesopt(rm)=mean(bayesopt_correct);
    %     stds_bayesopt(rm)=std(bayesopt_correct);
    end
    save(fullfile(S.path.prep,['response_model_predictions_' sname '.mat']), 'means_fitted','stds_fitted');
end

%totalmean = mean(means_fitted,1);
%totalstd = std(means_fitted,[],1);

% figure
% hold on
% bar(response_models,mean(means_fitted,1))
% errorbar(response_models,totalmean,totalstd,'.')
%plot(xlim,mean(means_bayesopt)*[1 1], 'k--')


switch S.resp_modelspec.responses{:}
    case 'Ch'
        % plot fraction predicted correctly
        figure
        hold on
        bar(response_models,mean(means_fitted,1))
        errorbar(response_models,mean(means_fitted,1),std(means_fitted,[],1),'.')
        %plot(xlim,mean(means_bayesopt)*[1 1], 'k--')

        % plot random corrected predictions
        figure
        hold on
        bar(response_models,mean(means_randcorr_fitted,1))
        errorbar(response_models,mean(means_randcorr_fitted,1),std(means_randcorr_fitted,[],1),'.')
    case 'RT'
        % plot correlation coefficients
        figure
        hold on
        bar(response_models,mean(cc,1))
        errorbar(response_models,mean(cc,1),std(cc,[],1),'.')
end
