%% Analysis: Perceptual model fitting

dbstop if error
clear all
close all
clear S

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'NLT'; %
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
S.select.subjects = {'e001'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {['Sequence_' S.expt '_OptionALPL_startblock*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
%S.select.blocks = {['Sequence_' S.expt '_OptionNLT_assoc*']}; % blocks to load (each a separate file) - empty means all of them, or not defined
%S.select.blocks = {['Sequence_' S.expt '_OptionALPL']}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.select.condtype = '';
S.path.datfile = ['C:\Data\MNP\participants\Participant_Data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add paths
run('C:\Data\Matlab\Matlab_files\MNP\MNP_addpaths');

% data import
S.load.prefixes = {'Output','Sequence'};
%[S,D]=SCIn_seq_import(S);
[S,D]=SCIn_data_import(S);

% unique save name extension
sname = datestr(now,30)
S.sname=sname;

% model fitting
S.prc_config = 'GBM_config'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_models=[3:8]; % 1 3 9 10 11 12
S.resp_model = [1]; 
for pm=1:length(S.perc_models)
    S.perc_model = S.perc_models(pm); 
    S=MNP_perceptual_models(S);
    S=MNP_response_models(S);
    
    % response data processing
    % NOW OCCURS AFTER MODEL SPECIFICATION TO ALLOW MULTIPLE ALPHAS
    S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
    S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
    S.accuracy.on = 1;
    S.RT.on = 1;
    S.RT.min = 0.5; % min RT to consider
    S.save.tables = 0;
    [S,D_prep]=SCIn_data_process(S,D);  
    
    % create HGF structure
    [D_prep,S] = MNP_prep4HGF(D_prep,S,1); % prepare response data for HGF
    
    % split data into training and testing sets (if we want to test for prediction of behaviour)
    S.frac_train = 0; % set to 0 to include all data in training set AND test set
    D_train=D_prep;
    if S.frac_train>0
        for d=1:length(D_prep)
            cond = D_prep(d).dt.design(2,:); % get conditions
            ucond = unique(cond);
            % random indices from each cond
            S.train_idx = [];
            for u = 1:length(ucond)
                cond_idx = find(cond==ucond(u));
                S.train_idx = [S.train_idx randsample(cond_idx,S.frac_train*length(cond_idx))];
            end
            S.test_idx = 1:length(D_train(d).HGF.u);
            S.test_idx(S.train_idx)=[];
            %D_train(d).HGF.u = D_prep(d).HGF.u(sort(S.train_idx),:);
            D_train(d).HGF.y(S.test_idx) = nan;
        end
    end

    S.HGF.plottraj = 0; % turn off if doing multiple simulations!
    D_fit=HGF_run(D_train,S,0);
    save(fullfile(S.path.hgf,'fitted',['MNP_fittedparameters_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']), 'D_fit', 'S');
    
end


% group model comparison
if 0
    close all
    S.perc_models=[1 3 9 10 11 12];
    S.resp_models = [1 2 3]; 
    S.sname='0_20180821T134505'; % 100% training
    %S.sname='0.5_20180818T070202'; % 50% training
    S.fname_pref= 'MNP_fittedparameters';
    S.fname_ext= ['_fractrain' S.sname '.mat'];
    [bmc.pm_out,bmc.pm_out]=HGF_group_model_comparison(S);
    
    % model freq and exceedence prob for families (both groups pooled)
    pm_ff  = bmc.pm_out{1,1}.families.Ef
    pm_fep = bmc.pm_out{1,1}.families.ep
    pm_ff  = bmc.pm_out{1,1}.families.Ef
    pm_fep = bmc.pm_out{1,1}.families.ep
end

% extract, tabulate and save parameters and summary stats of
% trajectories
if 0
    fname='MNP_fittedparameters_percmodel12_respmodel2_fractrain0_20180821T134505.mat';
    load(fname)
    [out.T,out.traj,out.param,out.rt] = MNP_extract_HGF_results(D_fit,S);
    [out.stats] = MNP_HGF_groupstatistics(out.T,{out.traj});
    save(fullfile(S.path.hgf,'fitted',strrep(fname,'fittedparameters','groupstatistics')), 'out');
    sname=strrep(strrep(fname,'fittedparameters','statstable'),'.mat','.xlsx');
    writetable(out.T,fullfile(S.path.hgf,'fitted',sname));
    
    % plot selected variable
    close all
    vs={'like_al0_1','like_al0_2','like_al0_3','like_al0_4','PL_om_2','PL_om_3','PL_sa_2_mean','PL_ud_1_mean','PL_ud_2_mean','PL_dau_mean','PL_epsi_2_mean','PL_epsi_3_mean'};
    [grps,~,grpind] = unique(out.T.groups);
    for v=1:length(vs)
        figure
        scatter(grpind,out.T.(vs{v}))
        set(gca,'XTick',1:2)
        set(gca,'XTickLabel',grps)
        title(vs{v})
    end
end

%MNP_sim_predict_behaviour(S,D_fit)