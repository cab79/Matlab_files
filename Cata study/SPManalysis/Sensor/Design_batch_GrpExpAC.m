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
D.batch_path = 'C:\Data\Matlab\Matlab_files\Cata study\SPManalysis\Sensor\Design_batch_GrpExpAC.m';
% template flexible factorial matlabbatch
D.ffbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_flexiblefactorial_template';
%  template SnPM matlabbatch
D.npbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_SnPM_template';
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
D.anapref = 't-3000_-2_b-3000_-2500'; %directory prefix for this specific analysis
D.subdirpref = '_mspm12_C'; % generic prefix for the SPM file type
%D.anapref = 'f10_t-3000_-2_b-3000_-2500'; %directory prefix for this specific analysis
%D.subdirpref = '_mrtf_spm12_C'; % generic prefix for the SPM file type
D.subdirsuff = '_orig_cleaned_SPNall'; % generic suffix for the EEGLAB analysis file
%D.subdirsuff = '_orig_cleaned_trialNmatch'; % generic suffix for the EEGLAB analysis file
D.folder =1; % Is the data in a subject-specific folder?
D.identifier=''; % optional identifer to add to end of outputted SPM folder name

% which codes to analyse in 'Include' columns in participant data file?
D.include_codes = [1];
% list of image names within each subject folder
D.imglist = {'scondition_c1.nii'
            'scondition_c3.nii'
            'scondition_c5.nii'
            'scondition_c7.nii'
            %'scondition_c2a.nii'
            %'scondition_c4a.nii'
            %'scondition_c6a.nii'
            %'scondition_c8a.nii'
            'scondition_c2b.nii'
            'scondition_c4b.nii'
            'scondition_c6b.nii'
            'scondition_c8b.nii'
            };
        
%% analysis design and parameters
%-------------------------------------------------------------
% specify parametric (1: spm) or non-parametric (2: SnPM) analysis. SnPM is
% limited to two-way interactions, so if three factors are select for
% interaction with SnPM it will output all 2-way interactions (actually, t-tests on subtracted data)
D.para = 1;
% specify a time window to analyse
%D.time_ana = [-3000 -500]; % applies a mask to the data
D.time_ana = [-2500 -1000]; % applies a mask to the data
% cond_list: each WITHIN SUBJECT factor (i.e. NOT including subject or group) is a column, each row is an
% image from imglist. Columns must be in same order as for 'factors' of type 'w' 
D.cond_list =  [
              1
              1
              1
              1
              2
              2
              2
              2
              %3
              %3
              %3
              %3
              ];
% factors and statistical model
D.factors = {'Grp', 'Exp', 'Subject'}; % must include a subject factor at the end; Group factor must be first if being used
D.factortype = {'g','w','s'}; % w = within, s = subject, g = subject group

% Main effects and interactions: 
%   - for spm, can specify the highest-level interaction to produc results
%   for all sub-interactions. Only main effects beyond those captured by
%   any interactions need to be listed, e.g. for Subject (only listed
%   Subject if there is no Group factor). E.g.
D.interactions = [1 1 0]; % one column per factor; one row per interaction
D.maineffects = [0 0 0]; % one column per factor 
%   - for snpm, only a single main effect or 2-way interaction can be performed each time, e.g.
%D.interactions = [0 0 0 0]; % one column per factor
%D.maineffects = [1 0 0 0]; % one column per factor 

D.grandmean = 0; % grand mean scaling value ('0' to turn off scaling)
D.globalnorm = 1; % Global normlisation: 1=off, 2 = proportional, 3 = ANCOVA

% names of nuisance covariates
%cov_names = {'Age','Gender'};
D.cov_names = {};

% the following are for spm analysis, not snpm
D.GMsca = [0 0 0]; %grand mean scaling
D.ancova = [0 0 0]; %covariate
% after model estimation, constrasts to display
D.fcontrasts = {
    [1 -1 -1 1], 'Grp * Exp'
    [1 1 -1 -1], 'Grp'
    [1 -1 1 -1], 'Exp'
    [0 -1 0 1], 'ExpC Grp'
    [-1 0 1 0], 'ExpA Grp'
    };

D.tcontrasts = {
    [0 1 0 -1], 'ExpC GrpA'
    [1 0 -1 0], 'ExpA GrpA'
    [0 -1 0 1], 'ExpC GrpB'
    [-1 0 1 0], 'ExpA GrpB'
    };
% the following are for SnPM, not SPM
D.nPerm = 5000; % permutations
D.vFWHM = [20 20 20]; % variance smoothing (should be same as data smoothing used)
D.bVolm = 1; % 1=high memory usage, but faster
D.ST_U = 0.001; % cluster forming threshold

% Write residuals? For normality tests
D.resid = 0;

%% run design_batch function
D=design_batch(D);

loadspm=0;
if D.para==1 && loadspm==1
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
