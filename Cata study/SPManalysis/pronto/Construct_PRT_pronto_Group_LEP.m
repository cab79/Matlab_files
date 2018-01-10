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
    D.grphead = {'Group'}; 
    D.inchead = 'Include';
% directory in which SPM analyses will be saved (new folder created)
D.spmstats_path = 'C:\Data\Catastrophising study\SPMstats\pronto';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix and suffix of subject folder names (within 'data_path') either side of subject ID
%D.anapref = 't-3000_-2_b-3000_-2500'; %directory prefix for this specific analysis
D.anapref = 't-5500_1500_b-5500_-5000'; %directory prefix for this specific analysis
D.subdirpref = '_mspm12_C'; % generic prefix for the SPM file type
%D.anapref = 'f10_t-3000_-2_b-3000_-2500'; %directory prefix for this specific analysis
%D.subdirpref = '_mrtf_spm12_C'; % generic prefix for the SPM file type
%D.subdirsuff = '_orig_cleaned_SPNall'; % generic suffix for the EEGLAB analysis file
D.subdirsuff = '_orig_cleaned_trialNmatch'; % generic suffix for the EEGLAB analysis file
D.folder =1; % Is the data in a subject-specific folder?
D.identifier='_GrpAvCond_gpc_ROI_perm100'; % optional identifer to add to end of outputted SPM folder name

% which codes to analyse in 'Include' columns in participant data file?
D.include_codes = [1];
% list of image names within each subject folder
D.imglist = {
            'scondition_c1.nii'
            'scondition_c3.nii'
            'scondition_c5.nii'
            'scondition_c7.nii'
            %'scondition_c2a.nii'
            %'scondition_c4a.nii'
            %'scondition_c6a.nii'
            %'scondition_c8a.nii'
            %'scondition_c2b.nii'
            %'scondition_c4b.nii'
            %'scondition_c6b.nii'
            %'scondition_c8b.nii'
            'scondition_c2.nii'
            'scondition_c4.nii'
            'scondition_c6.nii'
            'scondition_c8.nii'
            };
        
%% analysis design and parameters
%-------------------------------------------------------------
D.pronto = 1; % multivariate
% specify a time window to analyse
%D.time_ana = [-3000 -500]; % applies a mask to the data
%D.time_ana = [-2500 -1000]; % applies a temporal mask to the data (first level)
D.time_ana = [0 1500]; % applies a temporal mask to the data (first level)
D.timewin = [10];% apply windowing over the range of D.time_ana? Provide window size
% cond_list: each WITHIN SUBJECT factor (i.e. NOT including subject or group) is a column, each row is an
% image from imglist. Columns must be in same order as for 'factors' of type 'w' 
% For SnPM and Pronto, the second of two within-factors are subtracted.
D.cond_list =  [
              2
              2
              2
              2
              1
              1
              1
              1
              ];
D.grp_list = [1 2]; 
% factors and statistical model
D.factors = {'Grp','Exp','Subject'}; % must include a subject factor at the end; Group factor must be first if being used
D.factortype = {'g','w','s'}; % w = within, s = subject, g = subject group
%D.TrainTest = {[2 1],[],[]}; % select levels: model trains on first level and tests on second.
%D.grpcond = 1; % select a condition for the Grp Test

% Main effects and interactions: 
%   - for spm, can specify the highest-level interaction to produc results
%   for all sub-interactions. Only main effects beyond those captured by
%   any interactions need to be listed, e.g. for Subject (only listed
%   Subject if there is no Group factor). E.g.
D.interactions = [1 1 0]; % one column per factor; one row per interaction
D.maineffects = [0 0 0]; % one column per factor 
%   - for snpm and Pronto, only a single main effect or 2-way interaction can be performed each time, e.g.
%D.interactions = [0 0 0 0]; % one column per factor
%D.maineffects = [1 0 0 0]; % one column per factor 

% For SnPM or Pronto, if interactions enter the design then subtractions
% are needed by specifying 'contrast' here
D.fileoptype = 'meancond';
% overwrite previous images with same name
D.overwrite =0;

%D.grandmean = 0; % grand mean scaling value ('0' to turn off scaling)
%D.globalnorm = 1; % Global normlisation: 1=off, 2 = proportional, 3 = ANCOVA

% names of nuisance covariates - NOT IMPLEMENTED FOR PRONTO YET
%cov_names = {'Age','Gender'};
D.cov_names = {};

% the following are for spm analysis, not snpm
%D.GMsca = [0 0 0]; %grand mean scaling
%D.ancova = [0 0 0]; %covariate
% after model estimation, constrasts to display
%D.fcontrasts = {
%    };

% use kernel?
D.kernel = 1;

% permutation testing
D.permtest = 100;
D.saveallweights = 0; % Requires added code to prt_compute_weights_class.m, line 287:
                %if length(d.coeffs)>size(d.datamat,1)
                %    keep_idx = train_idx(find(ID(:,6)==1));
                %    d.coeffs = d.coeffs(train_idx);
                %end

% Data operations
% 1. Sample averaging (within blocks): constructs samples by computing the average of all
% scans within each block or event for each subject and condition.
% 2. Sample averaging (within subjects): constructs samples by computing the average of
% all scans within all blocks for each subject and condition.
% 3. Mean centre features using training data: subtract the voxel-wise mean from each
% data vector.
% 4. Divide data vectors by their norm: scales each data vector (i.e. each example) to lie
% on the unit hypersphere by dividing it by its Euclidean norm.
D.data_op  = {2 3}; 
% mean centring (same as 3 above)
%D.meancentre = 0;

% classification machine options:
%D.machine = 'svm_binary';
D.machine = 'gpc_binary';
%D.machine = 'gpc_multi'; % PRT_MODEL.M HAS BEEN MODIFIED TO ONLY ALLOW BINARY
%D.machine = 'mkl';

%D.cv_type = 'cv_lkso';D.nfolds = 2;
D.cv_type = 'cv_losgo';

%% run design_batch function
D=design_batch(D);

prt
