%function DCM_cata
% Reference: http://www.sciencedirect.com/science/article/pii/S105381191501037X#s0015
%--------------------------------------------------------------------------
clear all
close all
dbstop if error
%dbstop if error
spm('defaults','EEG');

% Use parallel pooling
S.parallel=1;

% Analyse small number of subjects/models in each group only? Speeds up script for debugging
% purposes
%S.run_subject_subsets = {1:8,9:16,17:24,25:34}; % subject subsets to analyse. Use '0' for all models
S.run_subject_subsets = {0}; % subject subsets to analyse. Use '0' for all subjects
S.run_model_subsets = 0; % model subsets to analyse. Use '0' for all models
S.run_num = 1; % label to assign to this series of runs
S.sub_ind = 0; % subject indices per group to analyse on this run
S.run_PEB = 1;
S.invert_all = 0; % invert all models (1), or only the full model using Bayesian Model Reduction for the other models (latter is faster - set to 0)
% prepare data or use the saved data (GCM struct)?
S.prepare_data = 1;
% save data from this run?
S.save_data = 0;
% exclude full model from BMC?
S.excl_full = 1; % automatically excludes if any A/B/C matrix values are <0
S.loadorsave=1;

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% root directory in which SPM data files are located
S.filepath = 'C:\Data\Catastrophising study\SPMdata'; 
% place to save source images
S.outpath = 'C:\Data\Catastrophising study\DCMdata\Timewin1a\Run1\modelA'; 
% name for group analysis output file
S.outname = 'DCM_GROUP'; % CURRENTLY UNUSED
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
S.pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
% names of headers in the above xls file:
    S.subhead = 'Subject';
    S.grphead = 'Group';
    S.inchead = 'Include';
%fiducials directory
S.fid_dir='C:\Data\Catastrophising study\meegfid';

%% specific file information for this analysis
%-------------------------------------------------------------
% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
S.fpref = 'mspm12';
S.fmid = ''; 
%S.fsuff = '_orig_cleaned.mat';
S.fsuff = '_orig_cleaned.mat';
% File containing location priors for dipoles
S.loc_file = 'loc.xlsx'; % xlsx file (in S.outpath) of locations with column 1: x, 2: y, 3: z, 4: location name
%model file name, in S.output
S.model_file = 'model.xlsx';

%% specific data settings for this analysis
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];
% groups to compare
S.grps = [1,2]; % use colon to analyse as separate groups or comma to analyse together (no group effect in model) 
% time and frequecy windows
%S.freqwin = []; % empty if not requiring freq analysis
S.timewin = [-5000 -4000]; 
S.basewin = [-5500 -5000]; 
% sort conditions into the order listed here prior to analysis
S.sortconds = {'c1','c2','c3','c4','c5','c6','c7','c8'};

% give a name to the within-subject contrast:
S.contrastname = 'None';
% which conditions to average together for this contrast?
S.conds = [1 1 1 1 1 1 1 1]; 
% Set the within-subject design: for condition effects 
% The number of trial-specific effects: one row per effect, ONE B MATRIX PER ROW
% ‘0 1’ to set the first trial type as the baseline for which to compare the second. 
% ‘0 1 1’ to set the first trial type as the baseline for which to compare both other conditions. 
% If there is no clear baseline condition, ‘-1 1’ to set the baseline as the average of two (MUST BE TWO CONDS ONLY TO USE THIS OPTION).
% To test a linear connectivity relationship, type all cond types in order
S.design = [0]; % for A matrix structure search
%S.design = [0 1]; % for B matrix condition effects

% between-subject covariate names in xls file
%S.cov_names = {'Age'};
S.cov_names = {};

%% Set up DCM model options
%==========================================================================

%--------------------------------------------------------------------------
% Parameters and options used for setting up model
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
%Model type: ’ERP’ is DCM for evoked responses; cross-spectral densities (CSD); induced responses (IND); phase coupling (PHA)
DCM.options.analysis = 'ERP'; % analyze evoked responses
%ERP is standard. CMC (default) and CMM are canonical microcircuit models based on the idea of predictive coding.
% reference: http://journal.frontiersin.org/article/10.3389/fncom.2013.00057/full
DCM.options.model    = 'CMC'; % ERP model
% Select single equivalent current dipole (ECD) for each source, or you use a patch on the 3D cortical surface (IMG).
DCM.options.spatial  = 'ECD'; % spatial model
% A projection of the data to a subspace is used to reduce the amount of data, using the 
% principal components of the prior covariance of the data
% Numebr of Modes = The number of principal components of the covariance matrix of the data. The default is 8.
% Should be greater than the number of expected sources.
DCM.options.Nmodes   = 13;     % nr of modes for data selection
%Select 0 for options.h to detrend and just model the mean. 
%Otherwise select the number of discrete cosine transform (DCT) terms you want to use to model low-frequency drifts (> 0). 
%DCT is similar to a fourier transform. 1 = 1 cycle of the cosine wave per time-period; 4 = 1 to 4 cycles.
DCM.options.h        = 0;     % nr of DCT components
% Onset parameter: to avoid attempting to model very early small deflections. 
% changing the onset prior has an effect on how your data are fitted. 
% default value of 60 ms onset time is a good value for many evoked responses where the first large deflection is seen around 100 ms. 
% However, this value is a prior, i.e., the inversion routine can adjust it. 
% Type several numbers (identical or not) only if there were several input stimuli (e.g. visual and auditory) – can be connected to different model sources.
DCM.options.onset    = -4900;    % selection of onset (prior mean)
% Duration (sd): makes it possible to vary the width of the input volley, separately for each of the inputs. 
% This can be used to model more closely the actual input structure (e.g. a long stimulus).
DCM.options.dur      = 100;    % Dispersion (sd)
DCM.options.D        = 4;     % downsampling - speeds up computation
DCM.options.han      = 0;     % Hanning removes the effect of beginning and end responses in time window.
% lock ECD orientations by introducing prior correlations. Useful for modelling bilateral symmetric sources (e.g., auditory cortices).
DCM.options.symmetry = 0;
% lock experimental effects by introducing prior correlations.
% ensures that all the changes in connectivity are the same. 
% This is useful when there is a specific hypothesis that some experimental factor increases (or decreases) all connection strengths.
DCM.options.lock = 0;
%“Optimise source locations” only works in combination with the “ECD” option and allows DCM more freedom with moving the dipoles as part of the optimisation process.
DCM.options.location = 0;
% Trial-specific inputs (C parameter estimation)
DCM.options.multiC = 0;

% maximum number of iterations [default = 64]
DCM.options.Nmax     = 32;

% 1: prepare data with forward model during inversion, 0: prepare the data beforehand
DCM.options.DATA     = 0;

run_DCM_group(S,DCM);
