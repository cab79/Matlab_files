% After entering the parameters, use icatb_batch_file_run(inputFile); 

dbstop if error
%% Modality. Options are fMRI and EEG
modalityType = 'EEG';

%% Type of analysis
% Options are 1, 2 and 3.
% 1 - Regular Group ICA
% 2 - Group ICA using icasso
% 3 - Group ICA using MST
which_analysis = 2;

%% ICASSO options.
% This variable will be used only when which_analysis variable is set to 2.
icasso_opts.sel_mode = 'randinit';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 10; % Number of times ICA will be run
% Most stable run estimate is based on these settings. 
icasso_opts.min_cluster_size = 2; % Minimum cluster size
icasso_opts.max_cluster_size = 10; % Max cluster size. Max is the no. of components

%% fMRI: Enter TR in seconds. If TRs vary across subjects, TR must be a row vector of length equal to the number of subjects.
%TR = [];

%% Group ica type
% Options are spatial or temporal for fMRI modality. Temporal for EEG. By default, spatial
% ica is run if not specified.
group_ica_type = 'temporal';

%% Parallel info
% enter mode serial or parallel. If parallel, enter number of
% sessions/workers to do job in parallel
parallel_info.mode = 'serial';
parallel_info.num_workers = 4;

%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m. 
% If you have selected option 3 (user specified settings) you need to manually set the PCA options. 
%
% Options are:
% 1 - Maximize Performance
% 2 - Less Memory Usage
% 3 - User Specified Settings
perfType = 1;

%% Design matrix selection
% Design matrix (SPM.mat) is used for sorting the components
% temporally (time courses) during display. Design matrix will not be used during the
% analysis stage except for SEMI-BLIND ICA.
% options are ('no', 'same_sub_same_sess', 'same_sub_diff_sess', 'diff_sub_diff_sess')
% 1. 'no' - means no design matrix.
% 2. 'same_sub_same_sess' - same design over subjects and sessions
% 3. 'same_sub_diff_sess' - same design matrix for subjects but different
% over sessions
% 4. 'diff_sub_diff_sess' - means one design matrix per subject.

keyword_designMatrix = 'no';

% specify location of design matrix here if you have selected 'same_sub_same_sess' or
% 'same_sub_diff_sess' option for keyword_designMatrix variable
OnedesignMat = '';

%% There are three ways to enter the subject data
% options are 1, 2, 3 or 4
dataSelectionMethod = 4;

%% Method 4
% Results will be saved in mat_results/runXXX/, where XXX is runIdx. Use a
% new runIdx for each dataset.
part='part2';
pth=['C:\Data\CORE\eeg\ana\prep\cleaned\' part];
files = dir(fullfile(pth,'*.set'));
for f = 1:length(files)
    erpFiles{f,1} = fullfile(pth,files(f).name);
end
baseLine = [-200 0];%Vector of length 2 to correct baseline
smoothTrials = 0;% - Option is provided to smooth trials using moving point average (moving_average_window).
moving_average_window = 0;%- Type of moving average window.
sorting_trials = 'No';%- Options are 'No', 'Condition & Latency'.
sort_idx = [];
dsample=[];

sort_idx = [
    5 6 13 14 21 22 ... %right hand, mismatch
    7 8 15 16 23 24 ... %right hand, standard
    1 2 9 10 17 18 ... %left hand, mismatch
    3 4 11 12 19 20]; %left hand, standard
outFileName = part;%- File name used when writing analyze or Nifti data.
outputDir = 'C:\Data\CORE\eeg\ana\groupICA\test';% - Output directory to save analyze or Nifti data.

% get trial indices
designfile = dir(fullfile('C:\Data\CORE\behaviour\raw',['dt_*' part '*']));
load(fullfile('C:\Data\CORE\behaviour\raw',designfile(1).name));
n_trials = size(dt.design,2); %max([mark(:).tnums]);

if ~exist(outputDir,'dir')
    mkdir(outputDir)
end
if 0
    CORE_icatb_eeg_import(char(erpFiles), baseLine, smoothTrials, moving_average_window, sorting_trials, sort_idx, outFileName, outputDir, dsample,n_trials)
end
% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of subjects
% and columns correspond to sessions. In the below example, there are 3
% subjects and 1 session. If you have multiple sessions, please see
% Input_data_subjects_2.m file.
input_data_file_patterns = strrep(erpFiles,'.set',['\' part '.mat']);
input_data_file_patterns = strrep(input_data_file_patterns,pth,outputDir);

% Input for design matrices will be used only if you have a design matrix
% for each subject i.e., if you have selected 'diff_sub_diff_sess' for
% variable keyword_designMatrix.
input_design_matrices = {};

% Enter no. of dummy scans to exclude from the group ICA analysis. If you have no dummy scans leave it as 0.
dummy_scans = 0;

%%%%%%%% End for Method 4 %%%%%%%%%%%%

%% Enter directory to put results of analysis
%outputDir = pth;

%% Enter Name (Prefix) Of Output Files
studyname = 'CORE';
prefix = [studyname '_' part];

%% Enter location (full file path) of the image file to use as mask
% or use Default mask which is []
maskFile = [];

%% Group PCA Type. Used for analysis on multiple subjects and sessions.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of timepoints between data-sets.
%
group_pca_type = 'subject specific';

%% Back reconstruction type. Options are str and gica
backReconType = 'gica';

%% Data Pre-processing options
% 1 - Remove mean per time point
% 2 - Remove mean per voxel
% 3 - Intensity normalization
% 4 - Variance normalization
preproc_type = 1;


%% PCA Type. Also see options associated with the selected pca option. EM
% PCA options and SVD PCA are commented.
% Options are 1, 2, 3, 4 and 5.
% 1 - Standard 
% 2 - Expectation Maximization
% 3 - SVD
% 4 - MPOWIT
% 5 - STP
pcaType = 1;

%% PCA options (Standard)

% a. Options are yes or no
% 1a. yes - Datasets are stacked. This option uses lot of memory depending
% on datasets, voxels and components.
% 2a. no - A pair of datasets are loaded at a time. This option uses least
% amount of memory and can run very slower if you have very large datasets.
pca_opts.stack_data = 'yes';

% b. Options are full or packed.
% 1b. full - Full storage of covariance matrix is stored in memory.
% 2b. packed - Lower triangular portion of covariance matrix is only stored in memory.
pca_opts.storage = 'full';

% c. Options are double or single.
% 1c. double - Double precision is used
% 2c. single - Floating point precision is used.
pca_opts.precision = 'double';

% d. Type of eigen solver. Options are selective or all
% 1d. selective - Selective eigen solver is used. If there are convergence
% issues, use option all.
% 2d. all - All eigen values are computed. This might run very slow if you
% are using packed storage. Use this only when selective option doesn't
% converge.

pca_opts.eig_solver = 'selective';


% %% PCA Options (Expectation Maximization)
% % a. Options are yes or no
% % 1a. yes - Datasets are stacked. This option uses lot of memory depending
% % on datasets, voxels and components.
% % 2a. no - A pair of datasets are loaded at a time. This option uses least
% % amount of memory and can run very slower if you have very large datasets.
% pca_opts.stack_data = 'yes';
% 
% % b. Options are double or single.
% % 1b. double - Double precision is used
% % 2b. single - Floating point precision is used.
% pca_opts.precision = 'single';
% 
% % c. Stopping tolerance 
% pca_opts.tolerance = 1e-4;
% 
% % d. Maximum no. of iterations
% pca_opts.max_iter = 1000;


% %% PCA Options (SVD)
% % a. Options are double or single.
% % 1a. double - Double precision is used
% % 2a. single - Floating point precision is used.
% pca_opts.precision = 'single';

% % b. Type of eigen solver. Options are selective or all
% % 1b. selective - svds function is used.
% % 2b. all - Economy size decomposition is used.
% pca_opts.solver = 'selective';


%% Maximum reduction steps you can select is 2. Options are 1 and 2. For temporal ica, only one data reduction step is
% used.
numReductionSteps = 2;

%% Batch Estimation. If 1 is specified then estimation of 
% the components takes place and the corresponding PC numbers are associated
% Options are 1 or 0
doEstimation = 0; 


%% MDL Estimation options. This variable will be used only if doEstimation is set to 1.
% Options are 'mean', 'median' and 'max' for each reduction step. The length of cell is equal to
% the no. of data reductions used.
estimation_opts.PC1 = 'mean';
estimation_opts.PC2 = 'mean';

%% Number of pc to reduce each subject down to at each reduction step
% The number of independent components the will be extracted is the same as 
% the number of principal components after the final data reduction step.  
numOfPC1 = min(92,icasso_opts.max_cluster_size+min(icasso_opts.max_cluster_size,10));
numOfPC2 = min(92,icasso_opts.max_cluster_size);

%% Scale the Results. Options are 0, 1, 2, 3 and 4
% 0 - Don't scale
% 1 - Scale to Percent signal change
% 2 - Scale to Z scores
% 3 - Normalize spatial maps using the maximum intensity value and multiply timecourses using the maximum intensity value
% 4 - Scale timecourses using the maximum intensity value and spatial maps using the standard deviation of timecourses
scaleType = 4;


%% 'Which ICA Algorithm Do You Want To Use';
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the
% command prompt.
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names

% 1 means infomax, 2 means fastICA, etc.
algoType = 1;


%% Specify atmost two reference function names if you select Semi-blind ICA algorithm.
% Reference function names can be acessed by loading SPM.mat in MATLAB and accessing 
% structure SPM.xX.name.
refFunNames = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};


%% Specify spatial reference files for constrained ICA (spatial) or gig-ica
refFiles = {which('ref_default_mode.nii'), which('ref_left_visuomotor.nii'), which('ref_right_visuomotor.nii')};

%% ICA Options - Name by value pairs in a cell array. Options will vary depending on the algorithm. See icatb_icaOptions for more details. Some options are shown below.
% Infomax -  {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0}
% FastICA - {'approach', 'symm', 'g', 'tanh', 'stabilization', 'on'}
icaOptions = {'posact', 'off', 'sphering', 'on', 'bias', 'on', 'extended', 0};

