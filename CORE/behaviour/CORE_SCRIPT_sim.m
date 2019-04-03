%% Analysis: Perceptual model fitting

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
S.select.groups = {};
S.select.subjects = {'CORE003'}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\CORE\participants\Participant_data_age_balanced.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

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
S.perc_model=[3]; % 1 3 9 10 11 12
S.resp_model = [4];
S=CORE_perceptual_models(S);
S=CORE_response_models(S);
S.HGF.plottraj = 1; % turn off if doing multiple simulations!

% which parameters?
%S.sim=[]; % specify here? or use generic parameters: S=CORE_sim_parameters(S)
%S=CORE_sim_parameters(S);
S.fitted_hgf = 'D_fit_r1_it13_n30.mat'; 
%S.fitted_hgf = 'CORE_fittedparameters_percmodel12_bayesopt_20181019T083824.mat';
condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3','PL_psi_1','PL_psi_2','PL_psi_3','PL_sa_1','PL_sa_2','PL_sa_3','PL_sahat_1','PL_sahat_2','PL_sahat_3','PL_mu_1','PL_mu_2','PL_mu_3','PL_muhat_1','PL_muhat_2','PL_muhat_3'};
conds = [1:24];
D_in.dt = D_prep.dt; 

S_in1=CORE_get_median_params_from_fitted(S,[1:15]); %group 1
D_sim(1)=HGF_run(D_prep,S_in1,1);
S_in1.condmean = condmean; S_in1.cond.condmean = conds;
D_in.HGF.fit = D_sim(1).HGF.sim;
[out.T(1,:),~,~,~] = CORE_extract_HGF_results(D_in,S_in1);

S_in2=CORE_get_median_params_from_fitted(S,[16:30]); %group 2
D_sim(2)=HGF_run(D_prep,S_in2,1);
S_in2.condmean = condmean; S_in2.cond.condmean = conds;
D_in.HGF.fit = D_sim(2).HGF.sim;
[out.T(2,:),~,~,~] = CORE_extract_HGF_results(D_in,S_in2);
%save(fullfile(S.path.hgf,'sim',['CORE_sim_percmodel' num2str(S.perc_model) '_' S.sname '.mat']), 'D_sim', 'S');

% replace parameters and simulate
S_sim = S_in2;
S_sim.sim.PL_om(2) = S_in1.sim.PL_om(2);
D_sim(3)=HGF_run(D_prep,S_sim,1);
S_sim.condmean = condmean; S_sim.cond.condmean = conds;
D_in.HGF.fit = D_sim(3).HGF.sim;
[out.T(3,:),~,~,~] = CORE_extract_HGF_results(D_in,S_sim);

S_sim = S_in2;
S_sim.sim.PL_om(3) = S_in1.sim.PL_om(3);
D_sim(4)=HGF_run(D_prep,S_sim,1);
S_sim.condmean = condmean; S_sim.cond.condmean = conds;
D_in.HGF.fit = D_sim(4).HGF.sim;
[out.T(4,:),~,~,~] = CORE_extract_HGF_results(D_in,S_sim);

S_sim = S_in2;
S_sim.sim.PL_om(2:3) = S_in1.sim.PL_om(2:3);
D_sim(5)=HGF_run(D_prep,S_sim,1);
S_sim.condmean = condmean; S_sim.cond.condmean = conds;
D_in.HGF.fit = D_sim(5).HGF.sim;
[out.T(5,:),~,~,~] = CORE_extract_HGF_results(D_in,S_sim);

% plot
y = [out.T.PL_epsi_1_condmean';
    out.T.PL_epsi_2_condmean';
    out.T.PL_epsi_3_condmean'];
x = categorical({'\epsilon_{1}','\epsilon_{2}','\epsilon_{3}'});

figure
cb = [1 0.25 0.25; 0.25 0.25 1; 255/255 204/255 0; 255/255 102/255 0; 0.75 0 0];
col = [1:5];
colormap(cb)
b=bar(x,y,'EdgeColor',[1 1 1]);
xtickangle(0)
for k = 1:size(y,2)
    set(b(k),'FaceColor',cb(col(k),:))
end
legend({'CRPS','HC','HC w/ CRPS \omega','HC w/ CRPS \theta','HC w/ CRPS \omega & \theta'})
box off
set(gca,'FontSize',18)