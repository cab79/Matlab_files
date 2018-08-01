% Analysis script for performing independent-samples tests over trials
% between two conditions.
% Input data can be epoched EEG .set files or .mat files with corresponding
% trial info files. Data can also be EEG or ICA components

clear all
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
restoredefaultpath

%% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

%% SET DATA PATHS/NAMES
clear S
S.path=struct;% clears the field
S.path.main = 'C:\Data\CORE\eeg\ana';
S.path.eeg = [S.path.main '\prep\cleaned\part2'];
S.path.stats = [S.path.main '\stats']; % folder to save outputs
S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; 
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

% layouts for cosmo
S.cfglayout = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.sfp';
S.cfgoutput = 'C:\Data\Matlab\Matlab_files\CORE\cosmo\modified functions\GSN92.lay';

% DATA TYPE (in order of speed)
% multicomp: all ICs considered as a group in a single multivariate analysis
% all_chan or recon: all ICs reconstructed into one set of channels; channels considered in a single multivariate analysis
% comp or chan: each IC or chan timecourse considered in separate analyses or as GFP
% comp_recon: each IC reconstructed into channels separately; each set of channels considered in a separate multivariate analyses
S.data_type='all_chan'; 

% ANALYSIS TYPE (see options at the top)
S.analysis_type='two_sample'; 

S.cond_idx = {
%     [1 2 9 10 17 18] %left hand, mismatch
%     [3 4 11 12 19 20] %left hand, standard
%     [5 6 13 14 21 22] %right hand, mismatch
%     [7 8 15 16 23 24]} %right hand, standard

    [1 2 9 10] %left hand, mismatch
    [3 4 11 12] %left hand, standard
    [5 6 13 14] %right hand, mismatch
    [7 8 15 16]}; %right hand, standard
S.contrast_rows = {[1 3],[2 4]}; % row of above cond_idx to contrast
S.total_samples = -200:799;
S.select_samples = 0:600;
S.smooth_samples = 10;
S.dsample = 4;
S.ndec=8; % trim data to a number of decimal places

% rank test
S.ranksum_on=0;

% TFCE settings (only use for recon and comp_recon, not on ICs)
S.tfce_on = 0; % analyses clusters over time on the GFP
S.cosmo_tfce_on = 0; % use cosmo if analysing clusters over channels/frequencies (+/- time)
S.tfce_tail = 1;
S.tfce_nperm=100; % for initial testing, between 100 and 1000. 5000 is most robust, but slow

% MVPA settings
S.mvpa_on=1;
S.SL_type = 'time';
S.search_radius = Inf;
S.use_measure = 'crossvalidation';
S.balance_dataset_and_partitions =1;
S.parti='nchunks'; % 'take-one-out', 'splithalf', 'oddeven', 'nchunks'
%S.use_classifier = 'LDA'; S.regularization = 0.5; S.matlab_lda = 0; S.logist = 0; S.output_weights = 1;
S.use_classifier = 'GP'; 
S.use_chunks = 'balance_targets';
S.nchunks=4;
S.average_train_count = 1;
S.average_train_resamplings = 1;

% savename
S.sname=datestr(now,30);

% process data
[S,D,stats] = CORE_eeg_trial_statistics(S);

% save stats
save(fullfile(S.path.stats,['stats_' S.analysis_type '_' S.sname '.mat']),'stats');

% plot data
if S.ranksum_on
    figure
    hold on
    bar(1:size(stats.ranksum.max,2),mean(stats.ranksum.max,1))
    errorbar(1:size(stats.ranksum.max,2),mean(stats.ranksum.max,1),std(stats.ranksum.max,[],1),'.')
    figure
    hold on
    bar(1:size(stats.ranksum.min,2),mean(stats.ranksum.min,1))
    errorbar(1:size(stats.ranksum.min,2),mean(stats.ranksum.min,1),std(stats.ranksum.min,[],1),'.')
end
if S.tfce_on
    for i=1:size(stats.tfce,3)
        figure
        hold on
        bar(1:size(stats.tfce,2),mean(stats.tfce(:,:,i),1))
        errorbar(1:size(stats.tfce,2),mean(stats.tfce(:,:,i),1),std(stats.tfce(:,:,i),[],1),'.')

        significant = find(stats.tfce(:,:,i)<=0.05)
        significant = S.select.subjects(significant)
    end
end
if S.mvpa_on
    figure
    imagesc(stats.mvpa_cv_acc,[0.4 0.6])
    colormap(parula(10))
    colorbar
    title([S.path.eeg ' ' S.analysis_type])
    figure
    hold on
    bar(1:size(stats.mvpa_cv_acc,2),mean(stats.mvpa_cv_acc,1))
    errorbar(1:size(stats.mvpa_cv_acc,2),mean(stats.mvpa_cv_acc,1),std(stats.mvpa_cv_acc,[],1),'.')
    plot(xlim,[0.6 0.6], 'k--')
    ylim([0.5 0.6])
    figure
    hold on
    bar(1:size(stats.mvpa_cv_acc,1),max(stats.mvpa_cv_acc,[],2))
    plot(xlim,[0.6 0.6], 'k--')
    ylim([0.5 0.7])
    
    disp(['mean = ' num2str(mean(stats.mvpa_cv_acc(:)))])
    disp(['std = ' num2str(std(stats.mvpa_cv_acc(:)))])
end

% group-level stats
if S.ranksum_on
    for c = 1:size(stats.signrank.stdgfp,1)
        stats.signrank.stdgfp(c,1)=signrank(stats.meandiff.stdgfp(:,c));
    end
end