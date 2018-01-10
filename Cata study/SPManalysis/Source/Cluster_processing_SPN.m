%% Cluster processing script 
% calls functions (need to be in Matlab search path): 
% Extract_clusters, Extract_cluster_waveforms, Convert_VOI_to_excel
clear all
dbstop if error
%% generic directories for all analyses for this study
%-------------------------------------------------------------
% load .xlsx file containing 'Subject', 'Group', and covariates
S.pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
% root directory in which source DATAFILES are located
S.datafile_path = 'C:\Data\Catastrophising study\SPMdata';
% root directory in which source IMAGES are located
S.data_path = 'C:\Data\Catastrophising study\SPMdata\sourceimages_SPNall_new';
% directory in which SPM analysis is saved 
S.spmstats_path = 'C:\Data\Catastrophising study\SPMstats\Source\1_grp\NoHanning\2nd_analysis_SPN';
% path to sensor space analysis where image_win.mat is saved
S.sensorpath = 'C:\Data\Catastrophising study\SPMstats\Include1\Between\t-3000_-2_b-3000_-2500_m_-2500_0_Grp_Exp_Subject_orig_cleaned_SPNall_spm';
% specific folder(s) containing the SPM stats for this analysis, 
% the original data file suffix,
% and the corresponding D.val (i.e. index of D.inv) from source analysis
S.spm_dir = {
    '_Time_Grp_Exp_Subject_spm_t-2346_-2266',1, '_orig_cleaned_SPNall.mat',{}
    %'_Time_Grp_Exp_Subject_spm_t-2500_-2332',1, '_orig_cleaned_SPNall.mat',{}
    %'_Time_Grp_Exp_Subject_spm_t-2296_-1950',1, '_orig_cleaned_SPNall.mat',{}
    %'_Time_Grp_Exp_Subject_spm_t-2144_-102',1, '_orig_cleaned_SPNall.mat',{}
    %'_Time_Grp_Exp_Subject_spm_t-1790_-1028',1, '_orig_cleaned_SPNall.mat',{}
    %'_Time_Grp_Exp_Subject_spm_t-938_-632',1, '_orig_cleaned_SPNall.mat',{}
    %'_Time_Grp_Exp_Subject_spm_t-350_-16',1, '_orig_cleaned_SPNall.mat',{}
};

%name of batch .mat file saved from design_batch.m and within same folder
%as SPM.mat
S.batch = 'matlabbatch.mat';
%name of subject information file in SPM directory (produced by
%design_batch)
S.subinfo = 'sub_info.mat';

%% for sourcewave extraction only:
%generic cluster image name for this run
%S.gclusname = 'comb_clus.nii';
%S.gclusname = 'VOI_c*.nii'; % Sensor
S.gclusname = 'c*_spm.nii'; % Source
% spm data file prefix for all analyses/images
S.prefix = 'mspm12_';
% using strsplit on the image filename, index of file part for subject
S.subname_index = 2;
% using strsplit on the image filename, index of file part for subject
%S.conname_index= 'end';
% Separate clusters for each unique value in each cluster image?
S.sep_clus =1;
% Use VOI data for MNI coordinates?
S.use_VOI = 0;
% scale outputs to VOI?
%s.scale_to_VOI = 1;
% min num of voxels required to constitutue a unique cluster region:
% smaller vallues to be deleted
S.clus_size_min = 100;
% min num of voxels required to constitutue a unique cluster region:
% smaller vallues to be joined to neighbouring regions
S.clus_size_min_join = 300;
% max number of non-continguous regions
%S.clus_obj_max = 1;
% analyse per subject/trial? (1: analyse all; 0: concatenate all)
S.ana_singlesub = 1;
% if SPM analysis has a 'time' factor (i.e. more than one time window),
% which level should be extracted?
S.timelev = 1;

%% for combining cluster data
S.use_aal = 1;
S.aal_path = 'C:\Data\Matlab\spm12\atlas\AAL2.nii';

%% for connectivity statistics only:
%  template SnPM matlabbatch
S.npbatch = 'C:\Data\Catastrophising study\SPMstats\matlabbatch_SnPM_template';
% directory in which image masks are saved
S.mask_path = 'C:\Data\Catastrophising study\SPMdata\masks';
% names of headers in the xls participant info file:
S.subhead = 'Subject';
S.grphead = 'Group';
S.inchead = 'Include';
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];

%% contrast, factor and level information
%-------------------------------------------------------------
% contrast name to process - must match that in Matlabbatch (i.e. from design-batch script)
% SUPERCEDED BY 4TH COLUMN OF SPM_DIR
S.contrasts={};
%S.contrasts={'Exp B Med'}; % leave empty to proccess ALL contrasts in Matlabbatch
%S.contrasts={'Int', 'Exp'}; % Post-stim
%S.contrasts={'T1 Grp * Exp','T1 Grp','T1 Exp'}; % Anticipation

S.tf =1; % 1 if F-contrast, 2 or T-contrast, blank if not using S.contrasts
% contrasts={'Exp'}; % example to enter one contrast only

% stats to save in a table for each contrast and cluster
S.clustab{1} = {'cluster','cluster','cluster','peak','peak','peak','','','',''; 
            'label','p(FWE-corr)','equivk','p(FWE-corr)','F','equivZ','x,y,z {mm}','x','y','z'}; 
        
S.clustab{2} = {'cluster','cluster','cluster','peak','peak','peak','','','',''; 
            'label','p(FWE-corr)','equivk','p(FWE-corr)','T','equivZ','x,y,z {mm}','x','y','z'}; 
        
% Factors and levels: for saving VOI information for later plotting
% 1: factor name in design matrix, 2: output factor name 3: factor levels. 
% Factors can be in order of desired output, not necessarily of
% input into SPM design matrix. Levels must be in the same order as SPM
% design, but characters don't need to match anything.
S.factlev = {
        %{'Att'},{'Attention Task'},{'No Task','Task'};
        {'Time'},{'LOI vs Baseline'},{'LOI','Baseline'};
        {'Grp'},{'Group'},{'High','Low'};
        %{'Att'},{'Attention'},{'Pain','Loc'};
        %{'Int'},{'Stimulus Intensity'},{'Low','Medium'};
        {'Exp'},{'Expectation Cues'},{'Low, Low','High, Low','High'};
        {'Subject'},{'Subject'},{}; % can leave Subject levels empty as these will be populated by sub_info file.
    };
S.subrow = 4; % row of above factlev containing the subject factor
S.timerow = 0; % row of above factlev containing the time (e.g. peak vs baseline) factor

% specific mask image (with fill path and extension) or leave empty
S.imgmask = '';
% cluster forming threshold
S.thresDesc = 'none'; % 'FWE' or 'none'
S.clusformthresh = 0.001;

%% setup the ROI network connectivity statistics
S.run_correlation_analysis = 0;
S.outputrank = 4; % set to 1000 to use highest possible rank / number of unique ROIs
S.rankminmax = 'max'; % use maximum or minimum rank
S.ranktolfactor = 1.1; % e.g. 1.1 increases tolerance by 10% each time
S.Regularize.do            = true;                           % use regularization on partial correlation matrices using the graphical lasso. 
S.Regularize.path          = 0.001;                          % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
S.Regularize.method        = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
S.Regularize.Prior.a = 1/3;
S.Regularize.Prior.b = 0;
S.Regularize.adaptivePath  = false;                          % adapt the regularization path if necessary
S.leakageCorrectionMethod  = 'symmetric';                      % choose from 'closest', 'symmetric', 'pairwise' or 'none'. 
S.nEmpiricalSamples        = 8;                              % convert correlations to standard normal z-statistics using a simulated empirical distribution. This controls how many times we simulate data of the same size as the analysis dataset
S.EnvelopeParams.windowLength = 1/5; % s                    % sliding window length for power envelope calculation. See Brookes 2011, 2012 and Luckhoo 2012. 
S.EnvelopeParams.useEnv = false; % envelope the data? CAB
S.EnvelopeParams.useFilter = true;                           % use a more sophisticated filter than a sliding window average
S.EnvelopeParams.takeLogs  = false;
S.frequencyBands           = {[0]};                           % a set of frequency bands for analysis. Set to empty to use broadband, 0 for ERP. The bandpass filtering is performed before orthogonalisation. 
S.paradigm = 'task';
%S.timecourseCreationMethod = 'PCA';                          % 'PCA', 'mean', 'peakVoxel' or 'spatialBasis'
%S.outputDirectory          = outDir;                         % Set a directory for the results output
S.statspackage = 'snpm'; % 'fsl' or 'snpm'
S.groupStatisticsMethod    = 'fixed-effects';                % 'mixed-effects' or 'fixed-effects' - only latter is an option with SnPM
%S.FDRalpha                 = 0.05;                           % false determination rate significance threshold
%S.sessionName              = {'sess1', 'sess2', 'sess3', 'sess4'}; 
%S.SaveCorrected.timeCourse    = false;
%S.SaveCorrected.envelopes      = false;
%S.SaveCorrected.variances     = false;
%S.SaveCorrected.ROIweightings = false;
S.SubjectLevel.conditionLabel   = {'Att1_Int1_Exp1', 'Att1_Int2_Exp1', 'Att1_Int1_Exp2', 'Att1_Int2_Exp2','Att2_Int1_Exp1', 'Att2_Int2_Exp1', 'Att2_Int1_Exp2', 'Att2_Int2_Exp2'};
S.SubjectLevel.designSummary    = {[1 0 0 0 0 0 0 0]', [0 1 0 0 0 0 0 0]', [0 0 1 0 0 0 0 0]', [0 0 0 1 0 0 0 0]',[0 0 0 0 1 0 0 0]', [0 0 0 0 0 1 0 0]', [0 0 0 0 0 0 1 0]', [0 0 0 0 0 0 0 1]'}; % summarise the design matrix by condition label
S.SubjectLevel.contrasts        = {[1 1 1 1 1 1 1 1]; [1 1 1 1 -1 -1 -1 -1]; [1 -1 1 -1 1 -1 1 -1];  [1 1 -1 -1 1 1 -1 -1]};  % each contrast is a new cell
S.SubjectLevel.interaction = [3 4]; % indicies of two contrasts to interact. First = paired test; second = subject-level contrast.
S.nSubs = 28; % reduce number of subjects to this number for paired tests if there is not enough memory
S.GroupLevel.designMatrix       = [ones(1,16), zeros(1,18); 
                                  zeros(1,16), ones(1,18)]';
S.SubjectLevel.subjectDesign    = eye(size(S.GroupLevel.designMatrix,1));
S.GroupLevel.contrasts          = [1  1;  % contrast 1
                                   1 -1]; % contrast 2

%% run functions (sit back and relax)

if isempty(strfind(S.spm_dir{1,1},'spm_t'))
    sd = S.spm_dir;
    load(fullfile(S.sensorpath,'image_win.mat'));
    for di = 1:size(image_win,1)-1 % first one is baseline
        tw = image_win{di+1,1};
        identifier = ['_t' num2str(tw(1)) '_' num2str(tw(2))];
        S.spm_dir{di,1} = [sd{1,1} identifier];
        S.spm_dir{di,2} = sd{1,2};
        S.spm_dir{di,3} = sd{1,3};
    end
end

Extract_clusters_source(S);
Convert_VOImat_to_excel(S);
%Extract_cluster_residuals(S);
%Normality_test_residuals(S);
%Combine_clusters_source(S);
%Extract_cluster_waveforms_source(S);
%SW_connectivity(S);
%SW_connectivity_results(S);