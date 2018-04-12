%% setup design batch for factorial design and estimation
clear all
%files required: participant data file with columns headed 'Subject' (should be characters, e.g. S1),
%'Group' (can be described by numbers or characters, but numbers recommended), 'Include' (must be numbers 
%with 0 meaning to exclude subject from analysis) and if you have covariates, columns headed with each
%covariate name. There should be just one header row, and one row for each
%subject.

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% name and location of the current design-batch file
D.batch_path = 'C:\Data\Matlab\Matlab_files\Cata study\SPManalysis\pronto\Construct_PRT_pronto.m';
% template flexible factorial matlabbatch
D.batch = 'C:\Data\Catastrophising study\SPMstats\pronto\PRT.mat';
% root directory in which subject-specific folders are located
D.data_path = 'C:\Data\Catastrophising study\SPMdata\sensorimages';
% directory in which image masks are saved
D.mask_path = 'C:\Data\Catastrophising study\SPMdata\masks';
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
D.pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
% names of headers in the above xls file:
    D.subhead = 'Subject';
    D.grphead = 'Group';
    D.inchead = 'Include';
% directory in which SPM analyses will be saved (new folder created)
D.spmstats_path = 'C:\Data\Catastrophising study\SPMstats\pronto';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix and suffix of subject folder names (within 'data_path') either side of subject ID
D.anapref = 't-3000_-2_b-3000_-2500'; %directory prefix for this specific analysis
D.subdirpref = '_mspm12_C'; % generic prefix for the SPM file type
%D.anapref = 'f10_t-3000_-2_b-3000_-2500'; %directory prefix for this specific analysis
%D.subdirpref = '_mrtf_spm12_C'; % generic prefix for the SPM file type
D.subdirsuff = '_orig_cleaned_SPNall'; % generic suffix for the EEGLAB analysis file
%D.subdirsuff = '_orig_cleaned_trialNmatch'; % generic suffix for the EEGLAB analysis file
D.folder =1; % Is the data in a subject-specific folder?
D.identifier='_gpr_MC_ROI'; % optional identifer to add to end of outputted SPM folder name

% which codes to analyse in 'Include' columns in participant data file?
D.include_codes = [1];
% list of image names within each subject folder
D.imglist = {
            'scondition_c1.nii'
            'scondition_c3.nii'
            'scondition_c5.nii'
            'scondition_c7.nii'
            'scondition_c2a.nii'
            'scondition_c4a.nii'
            'scondition_c6a.nii'
            'scondition_c8a.nii'
            'scondition_c2b.nii'
            'scondition_c4b.nii'
            'scondition_c6b.nii'
            'scondition_c8b.nii'
            };
        
%% analysis design and parameters
%-------------------------------------------------------------
D.pronto = 1; % multivariate
% specify a time window to analyse
%D.time_ana = [-3000 -500]; % applies a mask to the data
D.time_ana = [-2500 -2]; % applies a temporal mask to the data
D.timewin = [10];% apply windowing over the range of D.time_ana? Provide window size
% contrast: regression will take place on a single image per subject. So,
% if there are multiple images (multiple conditions) and new image will be
% calculated first based on a contrast. Use 1s and 2s, not -1 and 1.
D.regress_contrast =  [1 1 1 1 0 0 0 0 0 0 0 0];
%D.regress_contrast =  [1 1 1 1 1 1 1 1];
% unique name for this contrast to apply to images
%D.regress_contrastname =  'Mean';
D.fileoptype = 'meancond';

D.regress_contrastname =  'LowExp';
% which groups to include in regression? (empty = all subjects in one regression, otherwise names of groups to do separate regressions on)
D.regress_group = {};

% names of covariates in the regression
D.cov_names = {
    %'POMS1',
    %'Exp_effect',
    'Exp_effect_c0norm',
    };
% if there are many covariates, these will be included in one multiple
% regression, unless this next value is set to 1, in which case separate
% regressions will be run.
D.runmultiple = 1;
%D.cov_names = {};

% overwrite previous images with same name
D.overwrite =1;

% mean centring:
D.meancentre = 1;

% use kernel?
D.kernel = 1;

% permutation testing
D.permtest = 0;
D.saveallweights = 0;

% Data operations
% 1. Sample averaging (within blocks): constructs samples by computing the average of all
% scans within each block or event for each subject and condition.
% 2. Sample averaging (within subjects): constructs samples by computing the average of
% all scans within all blocks for each subject and condition.
% 3. Mean centre features using training data: subtract the voxel-wise mean from each
% data vector.
% 4. Divide data vectors by their norm: scales each data vector (i.e. each example) to lie
% on the unit hypersphere by dividing it by its Euclidean norm.
D.data_op  = {2}; 

% regression machine options:
%D.machine = 'svm_binary';
D.machine = 'gpr';
%D.machine = 'krr';


%% run design_batch function
D=design_batch(D);

prt
