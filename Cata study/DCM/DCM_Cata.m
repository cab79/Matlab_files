%function DCM_cata
% Reference: http://www.sciencedirect.com/science/article/pii/S105381191501037X#s0015
%--------------------------------------------------------------------------
clear all
close all
%dbstop if error
spm('defaults','EEG');

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

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% root directory in which SPM data files are located
S.filepath = 'C:\Data\Catastrophising study\SPMdata'; 
% place to save source images
S.outpath = 'C:\Data\Catastrophising study\DCMdata\_Time_Int_Exp_Subject_spm_t416_478\Int'; 
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
% File containing location priors for dipoles
S.loc_file = 'loc.xlsx'; % xlsx file (in S.outpath) of locations with column 1: x, 2: y, 3: z, 4: location name
%model file name, in S.output
S.model_file = 'modelA.xlsx';

%% specific data settings for this analysis
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];
% groups to compare
S.grps = [1,2]; % use colon to analyse as separate groups or comma to analyse together (no group effect in model) 
% time and frequecy windows
%S.freqwin = []; % empty if not requiring freq analysis
S.timewin = [0 800]; 
S.basewin = [-500 0]; 
% sort conditions into the order listed here prior to analysis
S.conds = {'c1','c2','c3','c4','c5','c6','c7','c8'};
% give a name to the within-subject contrast:
S.contrastname = 'Int';
% which conditions to contrast?
%S.contrast = [1 1 2 2 1 1 2 2]; 
S.contrast = [1 1 1 1 1 1 1 1]; 
% Set the "base condition":
% ‘0 1’ to set the first trial type as the baseline for which to compare the second. 
% ‘0 1 1’ to set the first trial type as the baseline for which to compare both other conditions. 
% If there is no clear baseline condition, ‘-1 1’ to set the baseline as the average of two (MUST BE TWO CONDS ONLY TO USE THIS OPTION).
% To test a linear connectivity relationship, type all cond types in order
S.basecond = [0]; 
%S.basecond = [0;1]; 
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
% reference: http://journal.frontiersin.org/article/10.3389/fncom.2013.00057/full
DCM.options.model    = 'ERP'; % ERP model
% Select single equivalent current dipole (ECD) for each source, or you use a patch on the 3D cortical surface (IMG).
DCM.options.spatial  = 'ECD'; % spatial model
% A projection of the data to a subspace is used to reduce the amount of data, using the 
% principal components of the prior covariance of the data
% Numebr of Modes = The number of principal components of the covariance matrix of the data. The default is 8.
% Should be greater than the number of expected sources.
DCM.options.Nmodes   = 14;     % nr of modes for data selection
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


%% RUN
cd(S.outpath)
%setup save name and run info
sii=1;
if exist(['Run_info' num2str(S.run_num) '.mat'],'file')
    load(['Run_info' num2str(S.run_num) '.mat']);
    for si = 1:length(S.run_subject_subsets)
        for su = 1:length(subs)
            ei(su) = isequal(subs(su),S.run_subject_subsets(si));
        end
        if ~any(ei)
            sii=si;
            break
        end
    end
    subs(1,length(subs)+1)=S.run_subject_subsets(sii);
else 
    subs=S.run_subject_subsets(sii);
end
save(['Run_info' num2str(S.run_num) '.mat'],'subs');
subjects = S.run_subject_subsets{sii}

if length(S.run_subject_subsets)>1
    S.sub_ind = S.run_subject_subsets{sii};
end

% options
DCM.options.Tdcm(1)  = S.timewin(1);     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = S.timewin(2);   % end of peri-stimulus time to be modelled
DCM.options.trials   = unique(S.contrast);
DCM.xU.X = S.basecond;
DCM.xU.name = {S.contrastname};
if size(S.basecond,2)>1
    S.basecond = S.basecond';
end

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------

% A matrix is the connections of the first trial/event
% NB: column index corresponds to the source area, and the row index to the
% target area.

% forward connections
DCM.A{1} = xlsread(S.model_file,'A_forward');
% backward connections
DCM.A{2} = xlsread(S.model_file,'A_backward');
% lateral connections
DCM.A{3} = xlsread(S.model_file,'A_lateral');
% B matrix: Gain modulations of connection strengths as set in the A-matrices
% Models the difference between the first and the other modelled evoked responses.
% For example, for two evoked responses, DCM explains the first response by
% using the A-matrix only. The 2nd response is modelled by modulating these connections 
% by the weights in the B-matrix.
% Set a single B model here to only test one model. 
%DCM.B{1} = DCM.A{1} + DCM.A{2};
%DCM.B{1}(1,1) = 1;
%DCM.B{1}(2,2) = 1;
DCM.B{1} = xlsread(S.model_file,'B');
% Alternatively, to test any combination of forward, backward, lateral and intrinsic connections in
% multipe models, type DCM.Bc = [1,2,3,4] for all, or for a subset. Leave
% empty to use the singel model defined above.
DCM.Ac = xlsread(S.model_file,'A_Acomb'); % order: F,B,L (no I)
DCM.Ac = find(DCM.Ac);
DCM.Anc = 1:3;
DCM.Anc(DCM.Ac) = [];
DCM.Bc = xlsread(S.model_file,'B_Acomb'); % order: F,B,L,I
DCM.Bc = find(DCM.Bc);
DCM.Bnc = 1:4;
DCM.Bnc(DCM.Bc) = [];
% C matrix: Inputs - can be to one or many areas
Cval = xlsread(S.model_file,'C');
DCM.C = Cval(1,:)';
% indices of nodes that should be tested for self-connection changes
DCM.Bi = Cval(2,:); 
DCM.Bi = find(DCM.Bi);

% read in locations data
[Lpos,Sname]=xlsread(fullfile(S.outpath,S.loc_file));
% specify the prior source locations (in mm in MNI coordinates). load the prior locations from a file
DCM.Lpos  = Lpos';
% enter the source names (one name in one row). 
DCM.Sname = Sname';
Nareas    = size(DCM.Lpos,2);

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

% create SubInd vector
SubInd_ana=cell2mat(SubInd);
% reduce subjects to those selected in this batch
if S.sub_ind
    SubInd_ana = SubInd_ana(S.sub_ind);
end

% between subject effects: constant, group difference
%--------------------------------------------------------------------------
if size(S.grps,1)==1
    Xb = [ones(sum(Nsub),1) ones(sum(Nsub),1)];
elseif size(S.grps,1)==2
    Xb = [ones(sum(Nsub),1) [ones(Nsub(1),1);-ones(Nsub(2),1)]];
end

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


B=[];
if any(DCM.Bc)
    % model space - within subject effects - creates B matrix of models to
    % compare
    % B matrix: Gain modulations of connection strengths as set in the A-matrices
    % In this case, we are creating 8 models: see Table 1 of http://www.sciencedirect.com/science/article/pii/S105381191501037X#s0015
    %----------------------------------------------------------------------
    % create (2^n x n) sparse matrix of indices permuted over n
    k = spm_perm_mtx(length(DCM.Bc));
    % for each reduced model, create sparse matrices of ones
    % sparse matrices are used when we expect most elements to contain zeros, and only a few non-zero elements.
    % this reduces memory size. In this case, it's done as an efficient way of
    % coding the connections (values of ones)
    for i = 1:(2^length(DCM.Bc));
        B{i}     = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
        % B comprises every combination (8 in total) of the options below: 
        try
            if k(i,1) && any(DCM.Bc==1)
                B{i} = B{i} + DCM.A{1}; % All forward connections
            end
        end
        try
            if k(i,2) && any(DCM.Bc==2)
                B{i} = B{i} + DCM.A{2}; % All backward connections
            end
        end
        try
            if k(i,3) && any(DCM.Bc==3)
                B{i} = B{i} + DCM.A{3}; % All lateral connections
            end
        end
        try
            if k(i,4) && any(DCM.Bc==4)
                B{i} = B{i} + sparse(DCM.Bi,DCM.Bi,1,Nareas,Nareas); % All intrinsic connections
            end
        end
        % add the zero elements in
        B{i}     = full(B{i});
    end
    nB = length(B);
elseif any(DCM.B{1})
    nB = sum(unique([DCM.B{:}])>0);
    B=[];
    k = spm_perm_mtx(nB);
    % for each possible combination of factors in B
    for i = 1:(2^nB);
        B{i}     = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
        % for each factor in B
        for f = 1:nB
            try
                if k(i,f)
                    B{i} = B{i} + double([DCM.B{:}]==f) + double([DCM.B{:}]==-1); 
                end
            end
        end
        if sum(B{i})==0
            B(i)=[];
        end
    end
else
    nB=1;
    B{1} = sparse(Nareas,Nareas);
end

%model space for A - must be no B matrix to use this
nA = sum(unique([DCM.A{:}])>0);
A=[];
if sum(DCM.B{:}(:))==0 && nA>1
    k = spm_perm_mtx(nA);
    % for each possible combination of factors in A
    for i = 1:(2^nA);
        %for each connection type in A
        for c = DCM.Ac
            A{i}{c} = sparse(Nareas,Nareas); % sparse 5x5 double array - for connecting 5 nodes
            % for each factor in A
            for f = 1:nA
                try
                    if k(i,f)
                        A{i}{c} = A{i}{c} + double(DCM.A{c}==f) + double(DCM.A{c}==-1); 
                    end
                end
            end
        end
        for c = DCM.Anc
            A{i}{c} = DCM.A{c}; 
        end
        if sum(A{i}{c})==0
            A(i)=[];
        end
    end
    nA = length(A);
else
    A{1}{1} = sparse(Nareas,Nareas);
    A{1}{2} = sparse(Nareas,Nareas);
    A{1}{3} = sparse(Nareas,Nareas);
    nA=1;
end

nC = sum(unique(DCM.C)>0);
C=[];
if nC>1
    k = spm_perm_mtx(nC);
    % for each possible combination of factors in A
    for i = 1:(2^nC);
        C{i}     = sparse(1,Nareas); % sparse 5x5 double array - for connecting 5 nodes
        % for each factor in C
        for f = 1:nC
            try
                if k(i,f)
                    C{i} = C{i} + double(DCM.C==f)'; 
                end
            end
        end
        if sum(C{i})==0
            C(i)=[];
        end
    end
    nC = length(C);
else
    C{1} = sparse(1,Nareas);
    nC=1;
end

% combinations of A, B and C
[Aind,Bind,Cind] = meshgrid(1:nA, 1:nB, 1:nC);
combs = [Aind(:),Bind(:),Cind(:)];
nMod = length(combs);

% loop through subjects to set up for inversion
%--------------------------------------------------------------------------
par=[];
if S.prepare_data==0 && exist('DCM_PEB_alldata.mat','file')
    load DCM_PEB_alldata
else
    subID={};
    files={};
    GCM=[];
    i=0;
    load(fullfile(S.fid_dir,'meegfid2.mat'));
    for i = 1:length(SubInd_ana)
        subID = pdata{SubInd_ana(i)+1,sub_col};
        if isnumeric(subID); subID = num2str(subID); end;
        fname = dir(fullfile(S.filepath,[S.fpref '*' subID '*' S.fmid  '*' S.fsuff]));
        files = fname.name;
        % Baseline Correction
        clear SB
        SB.D = fullfile(S.filepath,files);
        SB.timewin = S.basewin;
        SB.save = 0; % save in separate file
        SB.prefix = 'b'; % for output, only if saving a new file
        spm_eeg_bc(SB);

        % Data filename
        if isfield(DCM,'xY'); DCM = rmfield(DCM,'xY');end
        DCM.xY.Dfile = fullfile(S.filepath,files);

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

        for j = 1:nMod
            disp(['preparing data - subject number: ' subID ', model:' num2str(j)]); % display progress 
            if isfield(DCM,'M'); DCM = rmfield(DCM,'M');end
            DCM.name = ['DCM_' DCM.options.analysis '-ana_' DCM.options.model '-model_' DCM.options.spatial '-spatial_' num2str(j) '-model'];

            if j==1
                % Data and spatial model
                DCM  = spm_dcm_erp_data_CAB(DCM,1,S);
                xY =DCM.xY;
            else
                DCM.xY =xY;
            end

            DCM.A   = A{combs(j,1)};
            DCM.B   = B(combs(j,2));
            DCM.C   = C{combs(j,3)}';
            GCM{i,j} = DCM;
            
            % find number of free parameters to determine which is the full model
            [ii,rC,rE,Np] = spm_find_pC(GCM{i,j});
            par(:,j)     = sparse(ii,1,true,Np,1);
        end
        [~,fu] = max(sum(par));
        if i==1
            % Spatial model
            GCM{i,fu} = spm_dcm_erp_dipfit(GCM{i,fu},1);
            M.dipfit = GCM{i,fu}.M.dipfit;
        else
            GCM{i,fu}.M = M;
        end
    end
    if S.save_data
        disp('Saving data - can take a while...'); % display progress    
        save DCM_PEB_alldata GCM -v7.3
    end
end

% invert rreduced models (standard inversion) - takes HOURS!!!
%==========================================================================


%----------------------------------------------------------------------



if length(S.run_subject_subsets)>1
    sname = ['GCM_fit' num2str(S.run_num) num2str(S.run_subject_subsets{sii}(1)) '_' num2str(S.run_subject_subsets{sii}(end))];
else
    sname = ['GCM_fit' num2str(S.run_num)];
end
if S.invert_all
    if any(S.run_model_subsets>0)
        GCM_temp = GCM(:,S.run_model_subsets);
        GCM_temp = spm_dcm_fit(GCM_temp); % calls spm_dcm_erp for each subject and model
        GCM(:,S.run_model_subsets) = GCM_temp;
    else
        GCM = spm_dcm_fit(GCM); % calls spm_dcm_erp for each subject and model
    end
else
    %if length(S.run_subject_subsets)>1
    %    GCM_temp = spm_dcm_fit(GCM(S.run_subject_subsets{sii},fu)); % calls spm_dcm_erp for each subject and model
    %    GCM(S.run_subject_subsets{sii},fu) = GCM_temp;
    %else
        GCM(:,fu) = spm_dcm_fit(GCM(:,fu)); % calls spm_dcm_erp for each subject and model
    %end
end
subjects = S.run_subject_subsets{sii};
save([sname '.mat'],'GCM','DCM','Xb','subjects','A','B','C','combs','-v7.3');

if S.run_PEB
    lname = ['GCM_fit' num2str(S.run_num)];
    files = dir([lname '*.mat']);
    if length(files)>1
        if exist([lname '.mat'],'file')
            load([lname '.mat']);
            GCM_all=GCM;
            st=2;
        else
            st=1;
        end
        for f = st:length(files)
            load(files(f).name)
            GCM_all(subjects,:) = GCM;
            clear GCM
        end
        GCM=GCM_all;
        clear GCM_all
        save([sname '.mat'],'GCM','DCM','Xb','subjects','A','B','C','combs','-v7.3');
        for f = st:length(files)
            delete(files(f).name)
        end
    end
    loadorsave=1;
    DCM_PEB([lname '.mat'],S.outpath,S.run_num,fu,loadorsave)
end
