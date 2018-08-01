% Analysis script for Encoding and Decoding models
% This is a general script that can be tailored to specific analyses

% Encoding performs various forms of regression over trials (with or without cross-validation)
% The regression models for encoding can also be used to either generate corrected
% p-values (for cases where the tests are intended to be inferential) or to create cross-validated 
% parameter estimates.

% Encoding options:
% 1. Spearman correlation (SC) (with FDR and/or TFCE correction) - robust
% to outliers. TFCE only supports time (not chan) clustering presently.
% 2. linear multiple regression (MR) (with FDR and/or TFCE correction) -
% not robust to outliers, not robust to collinearity. TFCE only supports time (not chan) clustering presently.
% 3. non-hierarchical Bayesian linear regression: Parametric Empirical Bayes
% (PEB). This is used to obtain estimates of Log Model Evidence (LME) for model comparison, rather
% than for inferential tests against a null hypothesis.
% 4. Ridge regression (RR) for multiple collinear predictors. Includes
% cross-validation to select optimal lambda (penalisation). Only needed if
% there are a number of predictors in the model.
% 5. Bayesian regularised regression (BRR) with priors including Ridge and
% LASSO, useful for grouping together potentially collinear predictors,
% plus allows use of non-Gaussian distributions for robustness to outliers.
% Does not require cross-validation to optimise regularisation, but can
% test for predictive accuracy on test data if desired. Also outputs WAIC
% for model comparison.

% Decoding options:
% 1. BIEM: Bayesian inverted encoding model. Inverts the encoding model and cross-validates on test data. Currently implements the model published here:
    % Schoenmakers, S., Barth, M., Heskes, T., & van Gerven, M. (2013). Linear 
    % reconstruction of perceived images from human brain activity. NeuroImage,
    % 83, 951-961.
% 2. MVPA: Implements Gaussian Process Regression (GPR). 

% Input data can be epoched EEG .set files or .mat files with corresponding
% trial info files. Data can also be EEG or ICA components
% Predictor variables are RTs or HGF-estimated belief and prediction error
% trajectories.

clear all
close all
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
restoredefaultpath

%% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

%% SET DATA PATHS/NAMES: EEG .set files
clear S
S.path=struct;% clears the field
S.path.main = 'C:\Data\CORE\eeg\ana';
S.path.eeg = [S.path.main '\prep\cleaned\part2'];
S.path.stats = [S.path.main '\stats']; % folder to save outputs
S.path.hgf = ['C:\Data\CORE\behaviour\hgf\fitted\First draft\CORE_fittedparameters_percmodel6_respmodel3_20180710T131802.mat']; 
S.path.design = ['C:\Data\CORE\design']; % 
S.path.datfile = 'C:\Data\CORE\Participants\Participant_data.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
S.fname.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'set'}; 
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.load.suffixes = {'2_merged_cleaned'}; 
save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

%% SET DATA PATHS/NAMES: ICs
% clear S
% S.path=struct;% clears the field
% S.path.main = 'C:\Data\CORE\eeg\ana';
% S.path.eeg = [S.path.main '\groupICA\test'];
% S.path.stats = [S.path.main '\stats']; % folder to save outputs
% S.path.hgf = ['C:\Data\CORE\behaviour\hgf\fitted\First draft\CORE_fittedparameters_percmodel6_respmodel3_20180710T131802.mat']; 
% S.path.design = ['C:\Data\CORE\design']; % 
% S.path.datfile = 'C:\Data\CORE\Participants\Participant_data.xlsx'; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
% S.fname.parts = {'subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
% S.fname.ext = {'mat'}; 
% S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
% S.select.sessions = {};
% S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
% S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
% S.load.suffixes = {'2_merged_cleaned_grp-ica_c', '2_merged_cleaned_grp-fileinfo_'}; 
% save(fullfile(S.path.main,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable


% layouts for cosmo
S.cfglayout = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.sfp';
S.cfgoutput = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.lay';

% DATA TYPE (in order of speed)
% multicomp: all ICs considered as a group in a single multivariate analysis
% all_chan or recon: all ICs reconstructed into one set of channels; channels considered in a single multivariate analysis
% comp or chan: each IC or chan timecourse considered in separate analyses or as GFP
% comp_recon: each IC reconstructed into channels separately; each set of channels considered in a separate multivariate analyses
S.data_type='all_chan';

S.data_form = {'GFP','alldata'};

% ENCODING ANALYSIS TYPE (see options at the top) - leave empty for MVPA
S.analysis_type='RR'; 

S.cond_idx = {
    [1 2 9 10 17 18] %left hand, mismatch
    [3 4 11 12 19 20] %left hand, standard
    [5 6 13 14 21 22] %right hand, mismatch
    [7 8 15 16 23 24]} %right hand, standard

% which rows to subtract? (second cell must contain meaned data)
%S.row_subtract = {[1 3],[2 4]}; % mismatch trials minus mean of standards
%S.row_subtract = {}; % mismatch trials minus mean of standards

% rows of cond_idx to contrast
%S.contrast_rows = {}; % empty - all pooled into one (e.g. for regression)
%S.contrast_rows = {[1 3],[2 4]}; % e.g. fixed effects analysis
S.contrast_rows = {[1 3]}; % include mismatch only (e.g. if correlating with RT)
%S.contrast_rows = {[1:4]}; % include all (e.g. correlating with HGF traj)

% other data operations
S.total_samples = -200:799;
S.select_samples = 0:600;
S.smooth_samples = 10;
S.dsample = 4;
S.transform = 'arcsinh'; % arcsinh or notrans
S.zscore = 1;
S.ndec=8; % trim data to a number of decimal places
S.trainfrac = 0.5; % fraction of trials to use as training data (taken at random). Remainder used for testing. E.g. training for encoding vs. testing for decoding. Or, MPVA training vs. MVPA testing.
S.balance_conds =1; % random training trials balanced across conditions
S.num_runs = 1; % run analysis this many times, each time with a different random choice of training/test data

% Predictor variables
S.pred = 'RT'; % RT or HGF
% HGF trajectories, grouped
S.traj{1} = {
    {'PL'},{'mu','sa','muhat','sahat'};
    {'PR'},{'mu','sa','muhat','sahat'};
    }; % beliefs and their variance
S.traj{2} = {
    {'PL'},{'da','dau','ud','psi','epsi','wt'};
    {'PR'},{'da','dau','ud','psi','epsi','wt'};
    }; % prediction errors, updates and learning rates

% TFCE settings
S.tfce_on = 0; % for clustering over time only
S.cosmo_tfce_on = 0; % use cosmo if analysing clusters over channels/frequencies (+/- time). NB Cosmo does NOT implement regression/ccrrelation.
%S.tfce_test = S.analysis_type; % set to SC or MR above
S.tfce_tail = 2;
S.tfce_nperm=100; % for initial testing, between 100 and 1000. 5000 is most robust, but slow

% Ridge regression settings
S.rr.df_num = 100; % number of lambdas to do cross-val on
S.rr.folds = 5; % number of folds in the traindata for cross-validation
S.rr.z = 1; % determines which method to use. Either 1 or 0.
S.rr.sigma=0; % save sigma? RESULTS IN LARGE MATRICES AND MEMORY PROBLEMS

% Bayesian regularised regression (BRR) settings
S.brr.folds = 0;            % number of folds in the traindata. Set to 0 to not conduct predictive cross-validation.
S.brr.model = 'gaussian';   % error distribution - string, one of {'gaussian','laplace','t','binomial'}
S.brr.prior = 'ridge';        %- string, one of {'g','ridge','lasso','horseshoe','horseshoe+'}
S.brr.nsamples = 100;   %- number of posterior MCMC samples (Default: 1000)  
S.brr.burnin = 100;     %- number of burnin MCMC samples (Default: 1000)
S.brr.thin = 5;       %- level of thinning (Default: 5)
S.brr.catvars = [];    %- vector of variables (as column numbers of X) to treat
%                       as categorical variables, with appropriate expansion. 
%                       See examples\br_example5 (Default: none)
S.brr.nogrouping = false; %- stop automatic grouping of categorical predictors
%                       that is enabled with the 'catvars' options. (Default: false)
S.brr.usegroups = 0;     % ****Specified by S.traj cells**** - create groups of predictors. Grouping of variables
%                       works only with HS, HS+ and lasso prior
%                       distributions. The same variable can appear in
%                       multiple groups. See examples\br_example[9,10,11,12]  (Default: { [] } )  
S.brr.waic = true;       %- whether to calculate the WAIC -- disabling can lead
%                       to large speed-ups, especially for Gaussian models with large n
%                       (default: true)

% %BIEM settings
S.biem_on=1;
S.biem_prior = 'subject_training'; % options: 'group_training', 'subject_training', 'uniform'

% MVPA settings
S.mvpa_on=0;
S.SL_type = 'time';
S.search_radius = Inf; % data points, not ms.  Set to Inf to analyse all data points at once.
S.use_measure = 'crossvalidation';
S.balance_dataset_and_partitions =0; % turn on for classification, off for regression
S.parti='take-one-out'; % 'take-one-out', 'splithalf', 'oddeven', 'nchunks'
S.use_classifier = 'GP'; 
S.use_chunks = 'none'; % 'balance_targets' (over chunks) or 'none'
S.nchunks=0;
S.average_train_count = 1;
S.average_train_resamplings = 1;

% savename
S.sname=datestr(now,30);

% subject-level statistics (over trials)
[S,D,stats] = CORE_eeg_trial_statistics(S);

if 0
%plot SC
figure;imagesc(stats.spear.GFP(con).fdr_masked);
figure;imagesc(sum(stats.spear.alldata(con).fdr_masked,3));
% plot RR
imagesc(stats.RR.gfp.b{1, 1});
imagesc(squeeze(mean(stats.RR.alldata.b{1, 1},1)));

load('stats_MVPA_all_chan_RT_arcsinh_20180730T173749.mat'); mvpastats=stats;
load('stats_RR_all_chan_RT_arcsinh_20180730T174602.mat'); biemstats=stats;
ds = size(mvpastats.mvpa.alldata);

% plot MVPA predictions
figure; hold on
x=mean(reshape([mvpastats.mvpa.alldata.testdata_corr],ds(1),ds(2)),2);
y=mean(reshape([mvpastats.mvpa.alldata.samples],ds(1),ds(2)),2);
scatter(x,y); 
h=refline(1,0); h.Color = 'b';
line(xlim(), [0,0], 'Color', 'k');
line([0 0], ylim(), 'Color', 'k');
xlabel('test data predictive accuracy'); 
ylabel('training data predictive accuracy');

% plot MVPA predictions vs. BIEM
figure; hold on
x=mean(reshape([mvpastats.mvpa.alldata.testdata_corr],ds(1),ds(2)),2);
y=mean(reshape([biemstats.biem.alldata.rho],ds(1),ds(2)),2);
scatter(x,y); 
h=refline(1,0); h.Color = 'b';
line(xlim(), [0,0], 'Color', 'k');
line([0 0], ylim(), 'Color', 'k');
xlabel('MVPA test data predictive accuracy'); 
ylabel('BIEM test data predictive accuracy');

xticks=0:4:600;
load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat')
% plot MVPA weights
grp_weights=mean(cat(1,stats.mvpa(:).weights),1);
grp_weights = reshape(grp_weights,[],length(xticks));
[~,mi]=max(std(abs(grp_weights),[],1));
figure; imagesc(xticks,[],grp_weights); colormap jet; hold on; 
line(xticks([mi mi]),[1 92],'color','k','linewidth',2); title('weights')
figure; topoplot(grp_weights(:,mi),chanlocs); title('weights')

% plot MVPA transformed weights
grp_transweights=mean(cat(1,stats.mvpa(:).transweights),1);
grp_transweights = reshape(grp_transweights,[],length(xticks));
[~,mi]=max(std(grp_transweights,[],1));
figure; imagesc(xticks,[],grp_transweights); colormap jet; hold on; 
line(xticks([mi mi]),[1 92],'color','k','linewidth',2); title('transformed weights')
figure; topoplot(grp_transweights(:,mi),chanlocs); title('transformed weights')

% plot MVPA time-searchlight accuracy
xticks=0:4:600;
grp_acc=mean(cat(1,stats.mvpa(:).samples),1);
grp_std=std(cat(1,stats.mvpa(:).samples),[],1);
nsub = length(stats.mvpa); 
SEM = grp_std/sqrt(nsub);
tscore = -tinv(0.025,nsub-1);
CI = tscore*SEM;
upper = grp_acc+CI;
lower = grp_acc-CI;
figure; hold on
fill([xticks, fliplr(xticks)], [(upper), fliplr((lower))], ...
'b', 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
plot(xticks,grp_acc,'b'); 

end