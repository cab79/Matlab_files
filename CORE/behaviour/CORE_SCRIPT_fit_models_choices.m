%% Analysis: Perceptual model fitting and prediction of behaviour

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
[S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

% split data into training and testing sets
frac_train = 0;
if frac_train>0
    for d=1:length(D_prep)
        lengthtrain = round(frac_train*length(D_prep(d).HGF.u));
        D_train(d).HGF.u = D_prep(d).HGF.u(1:lengthtrain,:);
        D_test(d).HGF.u = D_prep(d).HGF.u(lengthtrain+1:end,:);
        D_train(d).HGF.y = D_prep(d).HGF.y(1:lengthtrain,:);
        D_test(d).HGF.y = D_prep(d).HGF.y(lengthtrain+1:end,:);
    end
else
    D_train=D_prep;
end

% model fitting
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_models=2:6;
S.sname = datestr(now,30);
for pm=1:length(S.perc_models)
    S.perc_model = S.perc_models(pm); 
    S=CORE_perceptual_models(S);
    S.resp_model = 1; 
    S=CORE_response_models(S);
    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    D_fit=HGF_run(D_train,S,0);
    save(fullfile(S.path.hgf,'fitted',['CORE_fittedparameters_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'D_fit');
    
    % add parameters to test data
    if frac_train>0
        for d=1:length(D_fit)
            D_test(d).HGF.fit.p_prc = D_fit(d).HGF.fit.p_prc;
        end
    else
        D_test=D_fit;
    end
    
    % Simulations
    S.numsimrep = 100; % number of simulations to run per parameter combination
    S.sim=[];
    [D_sim,S] = HGF_sim(D_test,S); 
    switch S.resp_modelspec.responses{:}
        case 'Ch'
            [means_fitted(:,pm),stds_fitted(:,pm),means_randcorr_fitted(:,pm),stds_randcorr_fitted(:,pm)]=HGF_test_model_predictions(D_test,D_sim,S);
            save(fullfile(S.path.hgf,'fitted',['CORE_fittedpredictions_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'means_fitted','stds_fitted');
        case 'RT'
            [cc(:,pm)]=HGF_test_model_predictions(D_test,D_sim,S);
            save(fullfile(S.path.hgf,'fitted',['CORE_fittedpredictions_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'cc');
    end
    
    % extract, tabulate and save parameters and summary stats of
    % trajectories
    [out.T,out.traj,out.param,out.rt] = CORE_extract_HGF_results(D_fit,S);
    [out.stats] = CORE_HGF_groupstatistics(out.T,{out.traj});
    save(fullfile(S.path.hgf,'fitted',['CORE_analysistables_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']), 'out');
end

% group model comparison
[out.bmc.posterior,out.bmc.out]=HGF_group_model_comparison(S)

switch S.resp_modelspec.responses{:}
    case 'Ch'
        % plot fraction predicted correctly
        figure
        hold on
        bar(S.perc_models,mean(means_fitted,1))
        errorbar(S.perc_models,mean(means_fitted,1),std(means_fitted,[],1),'.')
        %plot(xlim,mean(means_bayesopt)*[1 1], 'k--')

        % plot random corrected predictions
        figure
        hold on
        bar(S.perc_models,mean(means_randcorr_fitted,1))
        errorbar(S.perc_models,mean(means_randcorr_fitted,1),std(means_randcorr_fitted,[],1),'.')
    case 'RT'
        % plot correlation coefficients
        figure
        hold on
        bar(S.perc_models,mean(cc,1))
        errorbar(S.perc_models,mean(cc,1),std(cc,[],1),'.')
end

