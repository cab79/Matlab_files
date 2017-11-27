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
D.batch_path = 'C:\Data\Matlab\Matlab_files\Cata study\SPManalysis\Sensor\Design_batch_Regression.m';
% template flexible factorial matlabbatch
D.ffbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_regression_template';
%  template SnPM matlabbatch
D.npbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_SnPM_regression_template';
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
D.spmstats_path = 'C:\Data\Catastrophising study\SPMstats';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix and suffix of subject folder names (within 'data_path') either side of subject ID
D.anapref = 't-5500_-2500_b-5500_-5000'; %directory prefix for this specific analysis
%D.anapref = 't-3000_0_b-3000_-2500'; %directory prefix for this specific analysis
%D.anapref = 't-500_1500_b-500_0'; %directory prefix for this specific analysis
D.subdirpref = '_mspm12_C'; % generic prefix for the SPM file type
D.subdirsuff = '_orig_cleaned'; % generic suffix for the EEGLAB analysis file
%D.subdirsuff = '_orig_cleaned_trialNmatch'; % generic suffix for the EEGLAB analysis file
D.folder =1; % Is the data in a subject-specific folder?
D.identifier=''; % optional identifer to add to end of outputted SPM folder name

% which codes to analyse in 'Include' columns in participant data file?
D.include_codes = [1];
% list of image names within each subject folder
D.imglist = {'scondition_c1.nii'
            'scondition_c2.nii'
            'scondition_c3.nii'
            'scondition_c4.nii'
            'scondition_c5.nii'
            'scondition_c6.nii'
            'scondition_c7.nii'
            'scondition_c8.nii'};
        
%% analysis design and parameters
%-------------------------------------------------------------
% specify parametric (1: spm) or non-parametric (2: SnPM) analysis. SnPM is
% limited to two-way interactions, so if three factors are select for
% interaction with SnPM it will output all 2-way interactions (actually, t-tests on subtracted data)
D.para = 1;
% specify a time window to analyse
D.time_ana = []; % applies a mask to the data
% contrast: regression will take place on a single image per subject. So,
% if there are multiple images (multiple conditions) and new image will be
% calculated first based on a contrast. Use 1s and 2s, not -1 and 1.
D.regress_contrast =  [1 2 1 2 1 2 1 2];
%D.regress_contrast =  [1 1 1 1 1 1 1 1];
% unique name for this contrast to apply to images
%D.regress_contrastname =  'Mean';
D.regress_contrastname =  'Exp';
% Z-score normalise output image?
D.znorm = 1;
% which groups to include in regression? (empty = all subjects in one regression, otherwise names of groups to do separate regressions on)
D.regress_group = {};
D.fcontrasts = {
    [0 1], 'F'
    };

D.tcontrasts = {
    [0 1], 't+'
    [0 -1], 't-'
    };

% names of covariates in the regression
D.cov_names = {
    'Pain_Distraction','Pain_Reduction','Pain_Interfearnce_With_Task','pain_thresh',...
    'B1_upward_anxious','B1_upward_fear','B1_downward_anxious','B1_downward_fear',...
    'Task_Performance','Cue_Prediction','Cue_Attention','Cue_Influence','Prior_Expectancy_A',...
    'POMS1','Confidence_A','Exp_effect','Exp_effect_c0norm','Relief'
    };
% if there are many covariates, these will be included in one multiple
% regression, unless this next value is set to 1, in which case separate
% regressions will be run.
D.runmultiple = 1;
%D.cov_names = {};

D.grandmean = 0; % grand mean scaling value ('0' to turn off scaling)
D.globalnorm = 1; % Global normlisation: 1=off, 2 = proportional, 3 = ANCOVA

% the following are for spm analysis, not snpm
D.GMsca = [0 0 0 0]; %grand mean scaling
D.ancova = [0 0 0 0]; %covariate
% after model estimation, constrasts to display

% the following are for SnPM, not SPM
D.nPerm = 5000; % permutations
D.vFWHM = [20 20 20]; % variance smoothing (should be same as data smoothing used)
D.bVolm = 1; % 1=high memory usage, but faster
D.ST_U = 0.001; % cluster forming threshold

% Write residuals? For normality tests
D.resid = 0;

%% run design_batch function
if D.runmultiple
    D.cov_namesall=D.cov_names;
    for i = 1:length(D.cov_names)
        D.cov_names = D.cov_namesall(i);
        D=design_batch(D);
    end
else
    D=design_batch(D);
end

D.loadresults=0;
if D.para==1 && D.loadresults
    spm eeg
    load(fullfile(D.spm_path,'SPM.mat'));
    SPM.Im=[]; % no masking
    SPM.thresDesc = 'none'; % no correction
    SPM.u = 0.001; % uncorrected threshold
    SPM.k = 0; % extent threshold
    SPM.units = {'mm' 'mm' 'ms'};
    [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
    TabDat = spm_list('List',xSPM,hReg);
end
