%% Analysis: Perceptual model fitting
% This function approaches model inversion from an empirical Bayes
% perspective (based on VBA_MFX function in VBA toolbox), 
% whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% See: http://mbb-team.github.io/VBA-toolbox/wiki/VBA-MFX/

function CORE_condor_fit_step1

pth=pwd;

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'part2'; %
S.path.seq = fullfile(pth, 'Data'); % unprocessed data in original format
S.path.raw = fullfile(pth, 'Data'); % unprocessed data in original format
S.path.prep = fullfile(pth, 'Data'); % folder to save processed data
S.path.hgf = fullfile(pth, 'Data'); % folder to save processed data
S.path.design = fullfile(pth, 'Data'); % 
S.fname.parts = {'prefix','subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.groups = {};
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = fullfile(pth, 'Data','Participant_data.xlsx'); % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% unique save name extension
sname = datestr(now,30);
S.sname=sname;
if ~isempty(S.select.groups)
    S.sname = [S.sname '_' S.select.groups{:}];
end

% data import
S.load.prefixes = {'RT','dt'};
S.load.suffix = {'*'};
[S,D]=SCIn_data_import(S); 

% MFX options (see VBA-MFX)
opt=[];
opt = VBA_check_struct(opt, ...
    'TolFun'     , 2e-2  , ...     % Minimum change in the free energy
    'MaxIter'    , 32    , ...     % Maximum number of iterations
    'DisplayWin' , 1     , ...     % VB display window
    'verbose'    , 1       ...     % matlab window messages
    ) ;

% 1 - Model initialisation using pre-defined priors (ideally, Bayes optimal priors with wide variances)
S.prc_config = 'GBM_config_CORE_percmodel3_BOpriors'; S.obs_config = 'response_model_config'; S.nstim=[];S.bayes_opt=0;
S.perc_models=[3]; % 1 3 9 10 11 12
S.resp_models = [4]; % 2 12 21
S.perc_model = S.perc_models(1); 
S.resp_model = S.resp_models(1); 
S=CORE_perceptual_models(S);
S=CORE_response_models(S);

% MFX - Obtain and check priors from configs
% Default priors are used if priors are not explicitly provided through the
% priors_group structure. This means Gaussian(0,1) priors for the
% population mean of observation/evolution parameters and initial
% conditions, and Gamma(1,1) for the corresponding population precisions.
r.c_prc = eval([S.prc_config '(S,0)']); % runs the prc config function
r.c_obs = eval([S.obs_config '(r,S)']); % runs the obs config function
% convert log values to native
%[pvec.c_prc.priormus, pstruct] = GBM_transp(r, r.c_prc.priormus);
%[pvec.c_prc.priorsas, pstruct] = GBM_transp(r, r.c_prc.priorsas);
% put in VBA format for MFX analysis
priors_group.muPhi=r.c_obs.priormus';
priors_group.SigmaPhi=diag(r.c_obs.priorsas);
priors_group.b_vPhi=r.c_obs.priorsas';
priors_group.a_vPhi=double(r.c_obs.priorsas>0)'; % denominator needs to be either 1 (random variable) or inf (fixed effect)
priors_group.a_vPhi(priors_group.a_vPhi==0)=inf;
% put in VBA format for MFX analysis
priors_group.muTheta=r.c_prc.priormus';
priors_group.SigmaTheta=diag(r.c_prc.priorsas);
priors_group.b_vTheta=r.c_prc.priorsas';
priors_group.a_vTheta=double(r.c_prc.priorsas>0)'; % denominator needs to be either 1 (random variable) or inf (fixed effect)
priors_group.a_vTheta(priors_group.a_vTheta==0)=inf;
opt.priors_group = priors_group;

% obs: initialize the posterior on the population's mean and precision
p_group.muPhi = priors_group.muPhi;
p_group.SigmaPhi = priors_group.SigmaPhi;
iV_phi = VBA_inv(priors_group.SigmaPhi);
p_group.a_vPhi = priors_group.a_vPhi;
p_group.b_vPhi = priors_group.b_vPhi;
ind.phi_ffx = find(infLimit(p_group.a_vPhi,p_group.b_vPhi)==1);
ind.phi_in = find(diag(priors_group.SigmaPhi)~=0 & ~isnan(diag(priors_group.SigmaPhi)));
% prc: initialize the posterior on the population's mean and precision
p_group.muTheta = priors_group.muTheta;
p_group.SigmaTheta = priors_group.SigmaTheta;
temp=priors_group.SigmaTheta;
temp(isnan(temp)) = 0;
iV_theta = VBA_inv(temp); % uses temp because doesn't like nans
p_group.a_vTheta = priors_group.a_vTheta;
p_group.b_vTheta = priors_group.b_vTheta;
ind.theta_ffx = find(infLimit(p_group.a_vTheta,p_group.b_vTheta)==1);
ind.theta_in = find(diag(priors_group.SigmaTheta)~=0 & ~isnan(diag(priors_group.SigmaTheta)));

% response data processing
% NOW OCCURS AFTER MODEL SPECIFICATION TO ALLOW MULTIPLE ALPHAS
S.fitsim=1; % is the analysis of recorded data (1) or simulated data (2)?
S.meansim=0; % set to 1 to average results over repeated simulations (DONT USE FOR MODEL RECOVERY, ONLY FOR PLOTTING)
S.accuracy.on = 1;
S.RT.on = 1;
S.RT.min = 0.2; % min RT to consider
S.RT.max = 1;
S.RT.fill_non_response = 0; % RT value to add to trials in which there was no response: acts to include inaccurary info within RT data
S.save.tables = 0;
[S,D_prep]=CORE_data_process(S,D);  % specific function for CORE (bypasses SCIn_data_process)

% split data into training and testing sets (if we want to test for prediction of behaviour)
S.frac_train = 1; % set to 0 to include all data in training set AND test set
D_fit=D_prep;
if S.frac_train>0
    for d=1:length(D_prep)
        cond = D_prep(d).dt.design(2,:); % get conditions
        ucond = unique(cond);
        % random indices from each cond
        S.train_idx = [];
        for u = 1:length(ucond)
            cond_idx = find(cond==ucond(u));
            S.train_idx = [S.train_idx randsample(cond_idx,S.frac_train*length(cond_idx))];
        end
        S.test_idx = 1:length(D_fit(d).HGF.u);
        S.test_idx(S.train_idx)=[];
        %D_fit(d).HGF.u = D_prep(d).HGF.u(sort(S.train_idx),:);
        D_fit(d).HGF.y(S.test_idx) = nan;
    end
end

S.HGF.plottraj = 0; % turn off if doing multiple simulations!
S.parallel = 0; % parallelise over subjects

% invert
% if ~isempty(filebase)
%     if it_in>0
%         %it = it_in;
%         load(fullfile(S.path.hgf,filebase),'D_fit','mphi','Vphi','mtheta','Vtheta','p_group','F')
%     else
%         load(fullfile(S.path.hgf,filebase),'D_fit')
%     end
%     D_fit=D_fit(ismember({D_fit.subname},S.designmat(2:end,2)));
% else
    % obs: re-define within-subject priors
    r.c_obs.priormus=p_group.muPhi';
    %r.c_obs.priorsas=diag(p_group.SigmaPhi);
    r.c_obs.priorsas = transpose(p_group.b_vPhi./p_group.a_vPhi);
    S.c_obs = r.c_obs;
    % prc: re-define within-subject priors
    r.c_prc.priormus=p_group.muTheta';
    %r.c_prc.priorsas=diag(p_group.SigmaTheta);
    r.c_prc.priorsas = transpose(p_group.b_vTheta./p_group.a_vTheta);
    S.c_prc = r.c_prc; 

% end
S.priors_group = priors_group;
S.p_group = p_group;
S.ind = ind;
[nsub,nrand] = make_input_files(D_fit,S,pth); 

ns = sum(nsub);
create_job_submission_file(pth,ns,nrand)
    
% save info
fid = fopen(fullfile(pth,'nsub.txt'),'w');
fprintf(fid,'%d\n',nsub);
fclose(fid);
fid = fopen(fullfile(pth,'nrand.txt'),'w');
fprintf(fid,'%d',nrand);
fclose(fid);

quit 

function [nsub,nrand] = make_input_files(D,S,pth)

cd(pth)

prc_model = S.prc_config;
pmodel = strrep(prc_model,'_config',''); 
obs_model = S.obs_config;
opt_algo = 'tapas_quasinewton_optim_config';

S.use_y_col = find([any(strcmp(S.resp_modelspec.responses,'Ch')), any(strcmp(S.resp_modelspec.responses,'RT')),any(strcmp(S.resp_modelspec.responses,'EEG'))]);

if ~exist(fullfile(S.path.hgf,'input_files_orig'),'dir')
    mkdir(fullfile(S.path.hgf,'input_files_orig'))
end
if ~exist(fullfile(S.path.hgf,'input_files'),'dir')
    mkdir(fullfile(S.path.hgf,'input_files'))
end

for d = 1:length(D)
    
    disp(['creating input file, subject ' num2str(d) '/' num2str(length(D))])

    if ~isfield(S,'nstim')
        nst=length(D(d).HGF(1).u);
    else
        if isempty(S.nstim)
            nst=length(D(d).HGF(1).u);
        else
            nst=S.nstim;
        end
    end
    y = D(d).HGF(1).y(1:nst,:); 
    u = D(d).HGF(1).u(1:nst,:);
    
    save(fullfile(S.path.hgf,'input_files_orig',['input' num2str(d-1) '.mat']),'y','u','prc_model','obs_model','opt_algo','S')
        
end

% create randomised input files
load(fullfile(pth,'Data','grps_rand.mat'));
nrand = size(grps_rand,1);
%grp = importdata('grp.txt');
ngrp = length(unique(grps_rand(1,:)));
for g = 1:ngrp
    nsub(g) = sum(grps_rand(1,:)==g);
end
for i = 1:nrand
    disp(['creating randomisation ' num2str(i) '/' num2str(nrand)])
    ns=0;
    for g = 1:ngrp
        grp_ind = find(grps_rand(i,:)==g)-1;
        for gi = grp_ind
            ns = ns+1;

            % input file index
            ii = sum(nsub)*(i-1) + (ns-1);

            % copy file
            src=fullfile(S.path.hgf,'input_files_orig',['input' num2str(gi) '.mat']);
            des=fullfile(S.path.hgf,'input_files',['input' num2str(ii) '.mat']);
            copyfile(src,des);
        end
    end
end

disp('FINISHED');
    
function create_job_submission_file(pth,ns,nrand)

disp('creating job submission file')
A = {
    'executable=CORE_condor_fit_job.exe'
    'indexed_input_files=input.mat'
    'indexed_output_files=output.mat'
    'indexed_stdout=CORE_condor_fit_job.out'
    'indexed_stderr=CORE_condor_fit_job.err'
    'indexed_log=CORE_condor_fit_job.log'
    'max_run_time=60'
    ['total_jobs=' num2str(ns*nrand)]
};

fid = fopen(fullfile(pth, ['CORE_condor_fit_job_run.sub']),'w');
for i = 1:size(A,1)
    fprintf(fid,'%s\n',A{i});
end
fclose(fid);

function il = infLimit(a,b)
il = isinf(a).*eq(b,0);