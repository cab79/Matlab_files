%% Cluster processing script 
% calls functions (need to be in Matlab search path): 
% Extract_clusters, Extract_cluster_waveforms, Convert_VOI_to_excel
clear all
 
%% generic directories for all analyses for this study
%-------------------------------------------------------------
% load .xlsx file containing 'Subject', 'Group', and covariates
S.pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
% root directory in which subject-specific folders are located
S.data_path = 'C:\Data\Catastrophising study\SPMdata\sourceimages';
% directory in which SPM analysis is saved 
S.spmstats_path = 'C:\Data\Catastrophising study\SPMstats\Source';
% specific folder containing the SPM stats for this analysis
S.spm_dir = '_Time_Int_Exp_Subject_spm_t416_478';
%name of batch .mat file saved from design_batch.m and within same folder
%as SPM.mat
S.batch = 'matlabbatch.mat';
%name of subject information file in SPM directory (produced by
%design_batch)
S.subinfo = 'sub_info.mat';

%% contrast, factor and level information
%-------------------------------------------------------------
%contrast name to process - must match that in Matlabbatch (i.e. from design-batch script)
S.contrasts={'Int A'}; % leave empty to proccess ALL contrasts in Matlabbatch
S.tf =2; % 1 if F-contrast, 2 or T-contrast, blank if not using S.contrasts
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
        {'Time'},{'LOI vs Baseline'},{'LOI','Baseline'};
        {'Exp'},{'Expectation Cues'},{'Low, Low','High, Low'};
        {'Int'},{'Stimulus Intensity'},{'Low','Medium'};
        {'Subject'},{'Subject'},{}; % can leave Subject levels empty as these will be populated by sub_info file.
    };
S.subrow = 4; % row of above factlev containing the subject factor

% specific mask image (with fill path and extension) or leave empty
S.imgmask = '';
% cluster forming threshold
S.thresDesc = 'FWE'; % 'FWE' or 'none'
S.clusformthresh = 0.05;

%% run functions (sit back and relax)
Extract_clusters_source(S);
%Extract_cluster_waveforms(S);
Convert_VOImat_to_excel(S);
%Extract_cluster_residuals(S);
%Normality_test_residuals(S)