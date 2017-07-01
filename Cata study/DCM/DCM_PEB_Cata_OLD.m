%function DCM_PEB_cata
% Reference: http://www.sciencedirect.com/science/article/pii/S105381191501037X#s0015
%--------------------------------------------------------------------------
clear all
close all
dbstop if error
spm('defaults','EEG');

S.run_num = 1;
S.run_model_subsets = 2;

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% root directory in which SPM data files are located
S.filepath = 'C:\Data\Catastrophising study\SPMdata'; 
% place to save source images
S.outpath = 'C:\Data\Catastrophising study\DCMdata'; 
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
S.fsuff = '_orig_cleaned_trialNmatch.mat';

%% specific data settings for this analysis
% prepare data or use the saved data (GCM struct)?
S.prepare_data = 1;
% save data from this run?
S.save_data = 0;
% Analyse small number of subjects in each group only? Speeds up script for debugging
% purposes
S.sub_ind =1;
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];
% groups to compare
S.grps = [1 2]; 
% time and frequecy windows
%S.freqwin = []; % empty if not requiring freq analysis
S.timewin = [0 800]; 
S.basewin = [-500 0]; 
% sort conditions into the order listed here prior to analysis
S.conds = {'c1','c2','c3','c4','c5','c6','c7','c8'};
% give a name to the within-subject contrast:
S.contrastname = 'Int';
% which conditions to contrast?
S.contrast = [1 1 2 2 1 1 2 2]; 
% Set the "base condition":
% ‘0 1’ to set the first trial type as the baseline for which to compare the second. 
% ‘0 1 1’ to set the first trial type as the baseline for which to compare both other conditions. 
% If there is no clear baseline condition, ‘-1 1’ to set the baseline as the average of two (MUST BE TWO CONDS ONLY TO USE THIS OPTION).
% To test a linear connectivity relationship, type all cond types in order
S.basecond = [0;1]; 
% covariate names in xls file
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
DCM.options.model    = 'ERP'; % ERP model
% Select single equivalent current dipole (ECD) for each source, or you use a patch on the 3D cortical surface (IMG).
DCM.options.spatial  = 'IMG'; % spatial model
% A projection of the data to a subspace is used to reduce the amount of data, using use the 
% principal components of the prior covariance of the data; select the number of modes you wish to keep. The default is 8.
DCM.options.Nmodes   = 8;     % nr of modes for data selection
%Select 0 for options.h to detrend and just model the mean. 
%Otherwise select the number of discrete cosine transform (DCT) terms you want to use to model low-frequency drifts (> 0). 
%DCT is similar to a fourier transform. 1 = 1 cycle of the cosine wave per time-period; 4 = 1 to 4 cycles.
DCM.options.h        = 0;     % nr of DCT components
% Onset parameter: to avoid attempting to model very early small deflections. 
% changing the onset prior has an effect on how your data are fitted. 
% default value of 60 ms onset time is a good value for many evoked responses where the first large deflection is seen around 100 ms. 
% However, this value is a prior, i.e., the inversion routine can adjust it. 
% Type several numbers (identical or not) only if there were several input stimuli (e.g. visual and auditory) – can be connected to different model sources.
DCM.options.onset    = 300;    % selection of onset (prior mean)
% Duration (sd): makes it possible to vary the width of the input volley, separately for each of the inputs. 
% This can be used to model more closely the actual input structure (e.g. a long stimulus).
DCM.options.dur      = 150;    % Dispersion (sd)
DCM.options.D        = 4;     % downsampling - speeds up computation
DCM.options.han      = 1;     % Hanning removes the effect of beginning and end responses in time window.
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

%--------------------------------------------------------------------------
% Location priors for dipoles
%--------------------------------------------------------------------------
% specify the prior source locations (in mm in MNI coordinates). Can also load the prior locations from a file
DCM.Lpos  = [[-50; -10; 6] [50; 0; 12] [-34; -20; 62] [34; -12; 58] [-8; -48; 16]];
% enter the source names (one name in one row). 
DCM.Sname = {'left PI', 'right PI', 'left SP', 'right SP', 'PCC'};
Nareas    = size(DCM.Lpos,2);

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------

% A matrix is the connections of the first trial/event
% NB: column index corresponds to the source area, and the row index to the
% target area.

% forward connections
DCM.A{1} = zeros(Nareas,Nareas);
%DCM.A{1}(3,1) = 1;
%DCM.A{1}(4,2) = 1;
DCM.A{1}(5,1) = 1;
DCM.A{1}(5,2) = 1;
DCM.A{1}(5,3) = 1;
DCM.A{1}(5,4) = 1;

% backward connections
DCM.A{2} = zeros(Nareas,Nareas);
DCM.A{1}(1,5) = 1;
DCM.A{1}(2,5) = 1;
DCM.A{1}(3,5) = 1;
DCM.A{1}(4,5) = 1;

% lateral connections
DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(1,3) = 1;
DCM.A{3}(2,4) = 1;

% B matrix: Gain modulations of connection strengths as set in the A-matrices
% Models the difference between the first and the other modelled evoked responses.
% For example, for two evoked responses, DCM explains the first response by
% using the A-matrix only. The 2nd response is modelled by modulating these connections 
% by the weights in the B-matrix.
% Set a single B model here to only test one model. 
%DCM.B{1} = DCM.A{1} + DCM.A{2};
%DCM.B{1}(1,1) = 1;
%DCM.B{1}(2,2) = 1;
% Alternatively, to test any combination of forward, backward, lateral and intrinsic connections in
% multipe models, type DCM.Bw = [1,2,3,4] for all, or for a subset. Leave
% empty to use the singel model defined above.
DCM.Bw = [1]; % order: F,B,L,I
DCM.Bi = []; % name nodes that should be tested for self-connection changes


% C matrix: Inputs - can be to one or many areas
DCM.C = [1; 1; 0; 0; 0];


%% RUN

% options
DCM.options.Tdcm(1)  = S.timewin(1);     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = S.timewin(2);   % end of peri-stimulus time to be modelled
DCM.options.trials   = unique(S.contrast);
DCM.xU.X = S.basecond;
DCM.xU.name = {S.contrastname};

% NOTATION FOR SOME KEY FIELD NAMES:
%     DCM{i}.M.pE	- prior expectation of parameters
%     DCM{i}.M.pC	- prior covariances of parameters
%     DCM{i}.Ep		- posterior expectations
%     DCM{i}.Cp		- posterior covariance
%     DCM{i}.F		- free energy

if ~exist(S.outpath,'dir')
    mkdir(S.outpath);
end
cd(S.outpath)

[~,~,pdata] = xlsread(S.pdatfile);
grp_col = find(strcmp(pdata(1,:),S.grphead));
sub_col = find(strcmp(pdata(1,:),S.subhead));
inc_col = find(strcmp(pdata(1,:),S.inchead));

%Change char Grp inputs to numbers
grpdat = pdata(2:end,grp_col);
if isnumeric(grpdat{2,1})
    grptype = unique([grpdat{:}]);
    grptype(grptype==0)=[];
    Ngrp = length(grptype);
else
    grptype = unique(grpdat);
    grptype(isempty(grptype))=[];
    Ngrp = length(grptype);
    for g = 1:Ngrp
        grp_idx = cellfun(@(x) any(strcmp(grptype(g),x)), grpdat, 'UniformOutput', 0);
        grpdat(cell2mat(grp_idx)) = {[g]};
    end
end

% find index of subjects to include in analysis
SubInd = cell(Ngrp,1);
Subs = [];
inc_idx = cellfun(@(x) ismember(x,S.include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));

% find subject indices for each specific group
for g = 1:Ngrp
    grp_idx = find(cellfun(@(x) x==g, grpdat, 'UniformOutput', 1));
    SubInd{g,1} = intersect(inc_idx,grp_idx);
    Nsub(g,1) = length(SubInd{g,1});
end

% analyse group according to S.grps
SubInd_ana = {};
for g = 1:length(S.grps)
    SubInd_ana{g} = SubInd(S.grps(g));
end 
if S.sub_ind
    for g = 1:length(S.grps)
        SubInd_ana{g} = SubInd_ana{g}{1}(S.sub_ind);
    end 
end
Ngrp_ana = length(SubInd_ana);
Nsub_ana=[];
for g = 1:Ngrp_ana
    Nsub_ana(g) = length(SubInd_ana{g});
end

% between subject effects: constant, group difference
%--------------------------------------------------------------------------
Xb = [ones(sum(Nsub_ana),1) [ones(Nsub_ana(1),1);-ones(Nsub_ana(2),1)]];

% covariates
if ~isempty(S.cov_names)
    if ~strcmp(S.cov_names{1},'')
        for c = 1:length(S.cov_names)
            cov_col = find(strcmp(pdata(1,:),S.cov_names{c}));
            covdat = pdata(2:end,cov_col);

            % convert categorical (char) covariate data to numbers
            if ~isnumeric(covdat{2,1})
                covtype = unique(covdat);
                Ncovtype = length(covtype);
                for c = 1:Ncovtype
                    cov_idx = cellfun(@(x) any(strcmp(covtype(c),x)), covdat, 'UniformOutput', 0);
                    covdat(cell2mat(cov_idx)) = {[c]};
                end
            end
            covdat = covdat(vertcat(SubInd_ana{:}));

            Xb = horzcat(Xb,[covdat{:}]');
        end
    end
end

if ~isempty(DCM.Bw)
    % model space - within subject effects - creates B matrix of models to
    % compare
    % B matrix: Gain modulations of connection strengths as set in the A-matrices
    % In this case, we are creating 8 models: see Table 1 of http://www.sciencedirect.com/science/article/pii/S105381191501037X#s0015
    %----------------------------------------------------------------------
    % create (2^n x n) sparse matrix of indices permuted over n
    k = spm_perm_mtx(length(DCM.Bw));
    % for each reduced model, create sparse matrices of ones
    % sparse matrices are used when we expect most elements to contain zeros, and only a few non-zero elements.
    % this reduces memory size. In this case, it's done as an efficient way of
    % coding the connections (values of ones)
    for i = 1:(2^length(DCM.Bw));
        B{i}     = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
        % B comprises every combination (8 in total) of the options below: 
        try
            if k(i,1) && ~isempty(DCM.A{1})
                B{i} = B{i} + DCM.A{1}; % All forward connections
            end
        end
        try
            if k(i,2) && ~isempty(DCM.A{2})
                B{i} = B{i} + DCM.A{2}; % All backward connections
            end
        end
        try
            if k(i,3) && ~isempty(DCM.A{3})
                B{i} = B{i} + DCM.A{3}; % All lateral connections
            end
        end
        try
            if k(i,4)
                B{i} = B{i} + sparse(DCM.Bi,DCM.Bi,1,Nareas,Nareas); % All intrinsic connections
            end
        end
        % add the zero elements in
        B{i}     = full(B{i});
    end
else
    B=DCM.B;
end

% loop through subjects to set up for inversion
%--------------------------------------------------------------------------
if S.prepare_data==0 && exist('DCM_PEB_alldata.mat','file')
    load DCM_PEB_alldata
else
    subID={};
    files={};
    GCM=[];
    i=0;
    load(fullfile(S.fid_dir,'meegfid2.mat'));
    for g = 1:Ngrp_ana
        for s = 1:Nsub_ana(g)
            i=i+1;

            subID{g,s} = pdata{SubInd_ana{g}(s)+1,sub_col};
            if isnumeric(subID{g,s}); subID{g,s} = num2str(subID{g,s}); end;
            fname = dir(fullfile(S.filepath,[S.fpref '*' subID{g,s} '*' S.fmid  '*' S.fsuff]));
            files{g,s} = fname.name;
            % Baseline Correction
            clear SB
            SB.D = fullfile(S.filepath,files{g,s});
            SB.timewin = S.basewin;
            SB.save = 0; % save in separate file
            SB.prefix = 'b'; % for output, only if saving a new file
            spm_eeg_bc(SB);
            
            % Data filename
            if isfield(DCM,'xY'); DCM = rmfield(DCM,'xY');end
            DCM.xY.Dfile = fullfile(S.filepath,files{g,s});

            % prepare file: order conditions
            if isfield(S,'sortconds')
                if ~isempty(S.sortconds)
                    matlabbatch{1}.spm.meeg.preproc.prepare.D = {DCM.xY.Dfile};
                    matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = S.sortconds;
                    spm_jobman('initcfg')
                    spm_jobman('run',matlabbatch);
                end
            end

            % add fiducials
            load(DCM.xY.Dfile)
            if ~isequal(D.fiducials,meegfid)
                D.fiducials = meegfid;
                save(DCM.xY.Dfile,'D');
            end

            for j = 1:length(B)
                disp(['preparing data - subject number: ' num2str(i) ', model:' num2str(j)]); % display progress 
                if isfield(DCM,'M'); DCM = rmfield(DCM,'M');end
                DCM.B   = B(j);
                DCM.name = ['DCM_' DCM.options.analysis '-ana_' DCM.options.model '-model_' DCM.options.spatial '-spatial_' num2str(j) '-model'];

                if j==1
                    % Data and spatial model
                    DCM  = spm_dcm_erp_data(DCM);
                    xY =DCM.xY;
                else
                    DCM.xY =xY;
                end
                
                if i==1 && j==1
                    % Spatial model
                    DCM = spm_dcm_erp_dipfit(DCM,1);
                    M.dipfit = DCM.M.dipfit;
                else
                    DCM.M = M;
                end
                GCM{i,j} = DCM;
            end
        end

    end
    if S.save_data
        disp('Saving data - can take a while...'); % display progress    
        save DCM_PEB_alldata GCM -v7.3
    end
end

% invert rreduced models (standard inversion) - takes HOURS!!!
%==========================================================================
GCM_temp = GCM(:,S.run_model_subsets);
GCM_temp = spm_dcm_fit(GCM_temp); % calls spm_dcm_erp for each subject and model
load(['GCM_fit' num2str(S.run_num)]);
GCM(:,S.run_model_subsets) = GCM_temp;
save(['GCM_fit' num2str(S.run_num)],'GCM','DCM','GCM','Xb','Ns','-v7.3');

