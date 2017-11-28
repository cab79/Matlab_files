%% setup design batch function for factorial design and estimation
clear all
close all
%for fm = [2 3 17]
%files required: participant data file with columns headed 'Subject' (should be characters, e.g. S1),
%'Group' (can be described by numbers or characters, but numbers recommended), 'Include' (must be numbers 
%with 0 meaning to exclude subject from analysis) and if you have covariates, columns headed with each
%covariate name. There should be just one header row, and one row for each
%subject.

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% name and location of the current design-batch file
D.batch_path = 'C:\Data\Matlab\Matlab_files\Cata study\SPManalysis\Source\Design_batch_source_induced_SPN_GrpExp.m';
% template flexible factorial matlabbatch
D.ffbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_flexiblefactorial_template';
%  template SnPM matlabbatch
D.npbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_SnPM_template';
% root directory in which subject-specific folders are located
D.data_path = 'C:\Data\Catastrophising study\SPMdata\sourceimages_GS_1grp_noHan_SPN_induced';
% directory in which image masks are saved
D.mask_path = 'C:\Data\Catastrophising study\SPMdata\masks';
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
D.pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
% names of headers in the above xls file:
    D.subhead = 'Subject';
    D.grphead = 'Group';
    D.inchead = 'Include';
% directory in which SPM analyses will be saved (new folder created)
D.spmstats_path = 'C:\Data\Catastrophising study\SPMstats\Source\1_grp\NoHanning\2nd_analysis_SPN';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix and suffix of subject folder names (within 'data_path') either side of subject ID
D.anapref = ''; %directory prefix for this specific analysis
%D.anapref = 't-3000_0_b-3000_-2500'; %directory prefix for this specific analysis
%D.anapref = 't-500_1500_b-500_0'; %directory prefix for this specific analysis
D.subdirpref = 'spm12_C'; % generic prefix for the SPM file type
D.subdirsuff = '*'; % generic suffix for the EEGLAB analysis file
%D.subdirsuff = '_orig_cleaned'; % generic suffix for the EEGLAB analysis file
D.folder =0; % Is the data in a subject-specific folder?
% which codes to analyse in 'Include' columns in participant data file?
D.include_codes = [1];
% time identifer for source image files, to add to end of outputted SPM folder name
%identifier={'_t-2514_-2448','_t-2362_-2338','_t-2264_-2202','_t-2162_-2002','_t-1288_-1076'}; 
%D.baseline='_t-3000_-2500'; 
%identifier={'_t-4660_-4278','_t-4564_-4468','_t-4454_-4362','_t-4330_-4222','_t-3364_-3292','_t-3264_-3210','_t-3050_-2946','_t-2922_-2900','_t-2802_-2500'}; 
%D.baseline='_t-5500_-5000'; 
identifier={'_t-3000_0'}; 
D.baseline=''; 

for di = 1:length(identifier)
    D.identifier = identifier{di}

    %D.identifier = [D.identifier num2str(fm)];
    %D.data_path = [D.data_path num2str(fm)];

    % list of image names within each subject folder, or if not in a specific
    % folder then with a _ at the beginning
    D.imglist = {
                ['_*' D.identifier '_f8_12_c1.nii']
                ['_*' D.identifier '_f8_12_c2.nii']
                ['_*' D.identifier '_f8_12_c3.nii']
                ['_*' D.identifier '_f8_12_c4.nii']
                ['_*' D.identifier '_f8_12_c5.nii']
                ['_*' D.identifier '_f8_12_c6.nii']
                ['_*' D.identifier '_f8_12_c7.nii']
                ['_*' D.identifier '_f8_12_c8.nii']
                };

    %% analysis design and parameters
    %-------------------------------------------------------------
    % specify parametric (1: spm) or non-parametric (2: SnPM) analysis. SnPM is
    % limited to two-way interactions, so if three factors are select for
    % interaction with SnPM it will output all 2-way interactions (actually, t-tests on subtracted data)
    D.para = 1;
    % specify a time window to analyse
    D.time_ana = []; % applies a mask to the data
    % cond_list: each WITHIN SUBJECT factor (i.e. NOT including subject or group) is a column, each row is an
    % image from imglist. Columns must be in same order as for 'factors' of type 'w' 
    D.cond_list =  [
        1 1
        1 2
        1 1
        1 2
        2 1 
        2 2
        2 1 
        2 2
        ]';
    %D.cond_list =  [1 1 1 1]';
    % factors and statistical model
    D.factors = {'Att','Grp', 'Exp','Subject'}; % must include a subject factor at the end
    D.factortype = {'w','g','w','s'}; % w = within, s = subject, g = subject group

    % Main effects and interactions: 
    %   - for spm, can specify the highest-level interaction to produc results
    %   for all sub-interactions. Only main effects beyond those captured by
    %   any interactions need to be listed, e.g. for Subject (only listed
    %   Subject if there is no Group factor). E.g.
    D.interactions = [1 1 1 0]; % one column per factor; one row per interaction
    D.maineffects = [0 0 0 0]; % one column per factor 
    %   - for snpm, only a single main effect or 2-way interaction can be performed each time, e.g.
    %D.interactions = [0 0 0 0]; % one column per factor
    %D.maineffects = [0 0 1 0]; % one column per factor 

    % names of nuisance covariates
    %cov_names = {'Age','Gender'};
    D.cov_names = {};

    D.grandmean = 0; % grand mean scaling value ('0' to turn off scaling)
    D.globalnorm = 1; % Global normlisation: 1=off, 2 = proportional, 3 = ANCOVA

    % SPM only:
    D.GMsca = [0 0 0 0]; %grand mean scaling
    D.ancova = [0 0 0 0]; %covariate
    % after model estimation, constrasts to display (SPM, not SnPM)
     D.fcontrasts = {
        [1 -1 -1 1 -1 1 1 -1], 'Att * Grp * Exp'
        [-1 1 -1 1 1 -1 1 -1], 'Att * Exp'
        [1 1 -1 -1 -1 -1 1 1], 'Att * Grp'
        [-1 1 1 -1 -1 1 1 -1], 'Grp * Exp'
        [1 1 1 1 -1 -1 -1 -1], 'Att'
        [1 1 -1 -1 1 1 -1 -1], 'Grp'
        [-1 1 -1 1 -1 1 -1 1], 'Exp'
        };

    D.tcontrasts = {
        [1 1 1 1 -1 -1 -1 -1], 'Att A'
        [-1 -1 -1 -1 1 1 1 1], 'Att B'
        [1 1 -1 -1 1 1 -1 -1], 'Grp A'
        [-1 -1 1 1 -1 -1 1 1], 'Grp B'
        [-1 1 -1 1 -1 1 -1 1], 'Exp B'
        [1 -1 1 -1 1 -1 1 -1], 'Exp A'
        };

    % the following are for SnPM, not SPM
    D.nPerm = 5000; % permutations
    D.vFWHM = [20 20 20]; % variance smoothing (should be same as data smoothing used)
    D.bVolm = 1; % 1=high memory usage, but faster
    D.ST_U = 0.05; % cluster forming threshold

    % Write residuals? For normality tests
    D.resid = 0;

    %% run design_batch function
    D=design_batch(D);
end

%% load results
load_results=1;
if D.para==1 && load_results==1
    spm eeg
    load(fullfile(D.spm_path,'SPM.mat'));
    SPM.Im=[]; % no masking
    SPM.thresDesc = 'none'; % no correction
    SPM.u = 0.001; % uncorrected threshold
    SPM.k = 0; % extent threshold
    SPM.units = {'mm' 'mm' 'mm'};
    [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
    TabDat = spm_list('List',xSPM,hReg);
end
%end