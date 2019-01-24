%% Analysis: Perceptual model fitting
% This function approaches model inversion from an empirical Bayes
% perspective (based on VBA_MFX function in VBA toolbox), 
% whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% See: http://mbb-team.github.io/VBA-toolbox/wiki/VBA-MFX/

function CORE_SCRIPT_fit_response_models_empirical_priors(filebase,grp,it_in)

dbstop if error
close all

% FOLDER AND FILENAME DEFINITIONS
S.expt = 'part2'; %
S.path.seq = ['C:\Data\CORE\design']; % unprocessed data in original format
S.path.raw = ['C:\Data\CORE\behaviour\raw']; % unprocessed data in original format
S.path.prep = ['C:\Data\CORE\behaviour\processed']; % folder to save processed data
S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; % folder to save processed data
S.path.design = ['C:\Data\CORE\design']; % 
S.fname.parts = {'prefix','subject','suffix','ext'}; % parts of the input filename separated by underscores, e.g.: {'study','subject','session','block','cond'};
S.fname.ext = {'mat'}; 
S.select.groups = {grp};
S.select.subjects = {}; % either a single subject, or leave blank to process all subjects in folder
S.select.sessions = {};
S.select.blocks = {}; % blocks to load (each a separate file) - empty means all of them, or not defined
S.select.conds = {}; % conditions to load (each a separate file) - empty means all of them, or not defined
S.path.datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
%save(fullfile(S.path.prep,'S'),'S'); % saves 'S' - will be overwritten each time the script is run, so is just a temporary variable

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

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
S.parallel = 1; % parallelise over subjects

% 2 - evaluate free energies and parameter averages under the inital prior
% S.fname_pref= 'CORE_fittedparameters';
% S.fname_ext= ['_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat'];
% [LMEgrp,param_means]=HGF_group_model_comparison(S);
% from VBA_MFX
% F(1) = MFX_F(D_fit,p_group,priors_group,ind);

% for figures:
fig=figure('Name',grp);
dim.ns = length(D); % # subjects
dim.n_phi = length(p_group.muPhi);
dim.n_theta = length(p_group.muTheta);
dim.n = 0;
opt.dim = dim;
opt.priors_group = priors_group;
o_group.options = opt;
o_group.tStart  = tic;  % start time
p_sub = cell(dim.ns,1);
o_sub = cell(dim.ns,1);
o_group.initVBA.p_sub = p_sub;
o_group.initVBA.o_sub = o_sub;
[o_group.options] = VBA_displayMFX([],[],[],o_group,1,'off');


% 3 - iterate until convergence...
% We now update the within-subject effects as well as respective population
% moments by iteratively replacing the priors over within-subject effects by the VB
% estimate of the group mean and precision. The free energy of the ensuing
% procedure is computed for tracking algorithmic convergence.
stop = 0;
it = 1;
fprintf(1,['Iterative inversion...'])
while ~stop

    % invert
    if it==1 && ~isempty(filebase)
        if it_in>0
            it = it_in;
            load(fullfile(S.path.hgf,'fitted',filebase),'D_fit','mphi','Vphi','mtheta','Vtheta','p_group','F')
        else
            load(fullfile(S.path.hgf,'fitted',filebase),'D_fit')
        end
        D_fit=D_fit(ismember({D_fit.subname},S.designmat(2:end,2)));
        if strcmp(S.path.datfile,'C:\Data\CORE\participants\Participant_data_temp.xlsx')
            D_fit=D_fit(5:6);
        end
    else
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

        %fprintf(1,['Inital inversion with pre-defined priors...'])
        D_fit=HGF_run(D_fit,S,0); 
    end
    
    % store sufficient statistics
    optim_prc_ind = 1:length(ind.theta_in);
    optim_obs_ind = length(ind.theta_in)+1:length(ind.theta_in)+length(ind.phi_in);
    for d = 1:length(D_fit)
        mphi{it}(:,d) = D_fit(d).HGF.fit.p_obs.p;
        Vphi{it}{d} = nan(length(mphi{it}(:,d))); 
        Vphi{it}{d}(ind.phi_in,ind.phi_in) = D_fit(d).HGF.fit.optim.Sigma(optim_obs_ind,optim_obs_ind);
        % for mtheta, need to log values outputted from HGF
        [pvec, ~] = GBM_transp_log(r,D_fit(d).HGF.fit.p_prc.p);
        mtheta{it}(:,d) = pvec;
        % others are already ok as variances
        Vtheta{it}{d} = nan(length(mtheta{it}(:,d))); 
        Vtheta{it}{d}(ind.theta_in,ind.theta_in) = D_fit(d).HGF.fit.optim.Sigma(optim_prc_ind,optim_prc_ind);
    end
    
    % update moments of the parent population distribution
    [p_group.muPhi,p_group.SigmaPhi,p_group.a_vPhi,p_group.b_vPhi] = ...
            MFX_VBupdate(...
            priors_group.muPhi,...
            iV_phi,...
            mphi{it},...
            Vphi{it},...
            p_group.a_vPhi,...
            p_group.b_vPhi,...
            priors_group.a_vPhi,...
            priors_group.b_vPhi,...
            ind.phi_ffx,...
            ind.phi_in);
    [p_group.muTheta,p_group.SigmaTheta,p_group.a_vTheta,p_group.b_vTheta] = ...
            MFX_VBupdate(...
            priors_group.muTheta,...
            iV_theta,...
            mtheta{it},...
            Vtheta{it},...
            p_group.a_vTheta,...
            p_group.b_vTheta,...
            priors_group.a_vTheta,...
            priors_group.b_vTheta,...
            ind.theta_ffx,...
            ind.theta_in);
    
    % calculate free energy
    F(it) = MFX_F(D_fit,p_group,priors_group,ind)
    
    %figures
    o_group.F = F;
    o_group.it = it;
    o_group.ind = ind;
    [o_group.options] = VBA_displayMFX(p_sub,o_sub,p_group,o_group);
    
    %figures
%     figure(fig)
%     hold on
%     subplot(1,2,1);
%     scatter(it,F(it));
%     title('F')
%     drawnow
    
    % save
    save(fullfile(S.path.hgf,'fitted',['CORE_fittedparameters_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '_it' num2str(it) '.mat']),'D_fit','S','mphi','Vphi','mtheta','Vtheta','p_group','F');
    
    
    % stop rule
    if it>1
        dF = F(it) - F(it-1)
        %subplot(1,2,2);
        figure(fig)
        scatter(it,abs(dF));
        title('abs dF')
        %drawnow
        
        if abs(dF) <= opt.TolFun || it >= opt.MaxIter
            stop = 1;
        end
    end
    
    it = it +1;
    
end

% subfunctions (from VBA_MFX.m)
function [m,V,a,b] = MFX_VBupdate(m0,iV0,ms,Vs,a,b,a0,b0,indffx,indIn)
ns = size(ms,2);
n = size(m0,1);
sm = 0;
sv = 0;
wsm = 0;
sP = 0;
indrfx = setdiff(1:n,indffx);
indrfx = intersect(indrfx,indIn);
indffx = intersect(indffx,indIn);
iQ = diag(a(indrfx)./b(indrfx));
for i=1:ns
    % RFX
    sm = sm + ms(indrfx,i); % sum the means
    e = ms(indrfx,i)-m0(indrfx); % subject difference from prior mean
    sv = sv + e.^2 + diag(Vs{i}(indrfx,indrfx)); % add difference to variance
    % FFX
    tmp = VBA_inv(Vs{i});
    wsm = wsm + tmp*ms(:,i);
    sP = sP + tmp;
end
% RFX
V = zeros(n,n);
m = m0;
V(indrfx,indrfx) = VBA_inv(iV0(indrfx,indrfx)+ns*iQ); % updated SigmaTheta/Phi with prior variance
m(indrfx) = V(indrfx,indrfx)*(iV0(indrfx,indrfx)*m0(indrfx)+iQ*sm); % updated muTheta
a(indrfx) = a0(indrfx) + 0.5*ns; % resulting value always the same every iteration because a0 doesn't change
b(indrfx) = b0(indrfx) + 0.5*(sv+ns*diag(V(indrfx,indrfx))); % update prior variance given posteriors
% FFX
if ~isempty(indffx)
    tmp = VBA_inv(sP);
    V(indffx,indffx) = tmp(indffx,indffx);
    m(indffx) = V(indffx,indffx)*wsm(indffx);
end



function [F] = MFX_F(D,p_group,priors_group,ind)
% free energy computation
F = 0;
ns = length(D);
for i=1:ns
    F = F + D(ns).HGF.fit.optim.LME;
end
F = F + FreeEnergy_var(ns,...
    p_group.muPhi,p_group.SigmaPhi,...
    priors_group.muPhi,priors_group.SigmaPhi,...
    p_group.a_vPhi,p_group.b_vPhi,...
    priors_group.a_vPhi,priors_group.b_vPhi,...
    ind.phi_ffx,ind.phi_in);
F = F + FreeEnergy_var(ns,...
    p_group.muTheta,p_group.SigmaTheta,...
    priors_group.muTheta,priors_group.SigmaTheta,...
    p_group.a_vTheta,p_group.b_vTheta,...
    priors_group.a_vTheta,priors_group.b_vTheta,...
    ind.theta_ffx,ind.theta_in);

function F = FreeEnergy_var(ns,mu,V,mu0,V0,a,b,a0,b0,indffx,indIn)
% group-level variable-specific free energy correction term
n = length(mu);
indrfx = setdiff(1:n,indffx);
indrfx = intersect(indrfx,indIn);
n = length(indrfx);
e = mu(indrfx) - mu0(indrfx);
V = V(indrfx,indrfx);
V0 = V0(indrfx,indrfx);
a = a(indrfx);
b = b(indrfx);
a0 = a0(indrfx);
b0 = b0(indrfx);
iv0 = VBA_inv(V0);
F = -0.5*ns*sum(log(a./b)) ...
    + sum((a0+0.5*ns-1).*(psi(a)-log(b))) ...
    - sum((0.5*ns*diag(V)+b0).*a./b) ...
    + sum(a0.*log(b0) + gammaln(b0)) ...
    - 0.5*n*log(2*pi) ...
    - 0.5*VBA_logDet(V0) ...
    - 0.5*e'*iv0*e ...
    - 0.5*trace(iv0*V) ...
    + sum(entropyGamma(a,b)) + entropyGaussian(V) ...
    + 0.5*(ns-1).*length(indffx).*log(2*pi);

function S = entropyGamma(a,b)
S = a - log(b) + gammaln(a) + (1-a).*psi(a);

function S = entropyGaussian(V)
n = size(V,1);
S = 0.5*n*(1+log(2*pi)) + 0.5*VBA_logDet(V);

function il = infLimit(a,b)
il = isinf(a).*eq(b,0);

function combine_files
% combines two HGF results files from two separate group runs. Need to
% replace all the file names before running.
S.path.hgf = ['C:\Data\CORE\behaviour\hgf']; % folder to save processed data
Dbase = fullfile(S.path.hgf,'fitted');
Sfile = 'CORE_fittedparameters_percmodel3_respmodel4_fractrain0_20190109T180357.mat'; % file to get S.designmat (i.e. all subjects)
grp_files = {
    'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_CRPS_it14.mat'
    'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061503_HC_it14.mat'
    };
sname = 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it14.mat'; % save file name
load(fullfile(Dbase,Sfile),'S')
for g=1:length(grp_files)
    d=load(fullfile(Dbase,grp_files{g}),'D_fit');
    if g==1
        D_fit = d.D_fit;  
    else
        D_fit = [D_fit,d.D_fit];
    end
end
save(fullfile(Dbase,sname))