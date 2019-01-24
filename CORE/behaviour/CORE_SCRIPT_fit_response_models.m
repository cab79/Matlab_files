%% HGF model fitting script
% Runs multiple response models. Can run multiple times for different
% perceptual models. After all runs are complete, the end of the script
% allows factorial model comparison.

dbstop if error
clear all
close all

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
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
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

% model fitting
S.prc_config = 'GBM_config_CORE_percmodel2_BOpriors'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_models=[2]; 
S.resp_models = [7]; 
S.parallel = 1;
S.HGF.plottraj = 0; % turn off if doing multiple simulations!
for rm=1:length(S.resp_models)
    S.resp_model = S.resp_models(rm); 
    S.perc_model = S.perc_models(1); 
    S=CORE_perceptual_models(S);
    S=CORE_response_models(S);
    
    % response data processing
    % NOW OCCURS AFTER MODEL SPECIFICATION TO ALLOW MULTIPLE ALPHAS
    S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
    S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
    S.accuracy.on = 1;
    S.RT.on = 1;
    S.RT.min = 0.2; % min RT to consider
    S.RT.max = 1;
    S.RT.fill_non_response = 0; % RT value to add to trials in which there was no response: acts to include inaccurary info within RT data
    S.save.tables = 0;
    [S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

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

    D_fit=HGF_run(D_train,S,0);
    save(fullfile(S.path.hgf,'fitted',['CORE_fittedparameters_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']), 'D_fit', 'S');
    
end


% group model comparison
if 0
    close all
    S.sname='0_20190112T081048'; 
    S.perc_models=[1 2 3];
    S.resp_models = [1 2 3 4 5 6 7]; 
    %S.perc_models=[1205 1210 1215 1220 1225 1225 1230 1235 1240 1245 1250]; S.sname='0_20190101T162309'; % alpha
    %S.perc_models=[121 122 124 128 1216 1232 1264]; S.sname='0_20181228T104344'; % alpha
    %S.perc_models=[1205 121 122]; S.sname='0_20181221T190047'; % om2
    %S.perc_models=[121]; S.sname='0_20181227T172306'; 
    %S.resp_models = [15]; 
    
    S.fname_pref= 'CORE_fittedparameters';
    S.fname_ext= ['_fractrain' S.sname '.mat'];
    S.pep_flag = 1; % display PEP instead of EP
    S.family_on = 1;
    [~,~,bmc.pm_out,bmc.rm_out]=HGF_group_model_comparison(S);
    
    % model freq and exceedence prob for families (both groups pooled)
    pm_ff  = bmc.pm_out{1,1}.families.Ef
    pm_fep = bmc.pm_out{1,1}.families.ep
    rm_ff  = bmc.rm_out{1,1}.families.Ef
    rm_fep = bmc.rm_out{1,1}.families.ep
end
