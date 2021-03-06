%% Analysis: Perceptual model fitting
% This function approaches model inversion from an empirical Bayes
% perspective (based on VBA_MFX function in VBA toolbox), 
% whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% See: http://mbb-team.github.io/VBA-toolbox/wiki/VBA-MFX/

function CORE_condor_fit_step2

pth = pwd;

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

% get input and output file info
[in,out] = CORE_condor_monitor_outputs;

nsub = importdata('nsub.txt');
nrand = importdata('nrand.txt');
it = importdata('it.txt');

for nr = 1:nrand
    D_fit=struct;
    
    for d = 1:nsub
        
        % input file index
        ii = nrand*(nr-1) + (d-1);
        
        % create D_fit
        in_fname = ['input' num2str(ii) '.mat'];
        out_fname = ['output' num2str(ii) '.mat'];
        load(in_fname,'u','y');
        fit=load(out_fname,'out');
        D_fit(d).HGF.u=u;
        D_fit(d).HGF.y=y;
        D_fit(d).HGF.fit=fit.out;
    end

    % get group prior stats
    ii = nrand*(nr-1);
    load(['input' num2str(ii) '.mat'],'S');
    priors_group = S.priors_group;
    p_group = S.p_group;
    ind = S.ind;
    iV_phi = VBA_inv(priors_group.SigmaPhi);
    temp=priors_group.SigmaTheta;
    temp(isnan(temp)) = 0;
    iV_theta = VBA_inv(temp); % uses temp because doesn't like nans

    % save D_fit,S
    save_dir = fullfile(pwd,'Data/fitted');
    if ~exist(save_dir)
        mkdir(save_dir)
    end
    save(fullfile(save_dir,['D_fit_r' num2str(nr) '_it' num2str(it)]),'D_fit','S')

    % store sufficient statistics
    r.c_prc = eval([S.prc_config '(S,0)']); % runs the prc config function
    r.c_obs = eval([S.obs_config '(r,S)']); % runs the obs config function
    optim_prc_ind = 1:length(ind.theta_in);
    optim_obs_ind = length(ind.theta_in)+1:length(ind.theta_in)+length(ind.phi_in);
    for d = 1:length(D_fit)
        mphi(:,d) = D_fit(d).HGF.fit.p_obs.p;
        Vphi{d} = nan(length(mphi(:,d))); 
        Vphi{d}(ind.phi_in,ind.phi_in) = D_fit(d).HGF.fit.optim.Sigma(optim_obs_ind,optim_obs_ind);
        % for mtheta, need to log values outputted from HGF
        [pvec, ~] = GBM_transp_log(r,D_fit(d).HGF.fit.p_prc.p);
        mtheta(:,d) = pvec;
        % others are already ok as variances
        Vtheta{d} = nan(length(mtheta(:,d))); 
        Vtheta{d}(ind.theta_in,ind.theta_in) = D_fit(d).HGF.fit.optim.Sigma(optim_prc_ind,optim_prc_ind);
    end

    % update moments of the parent population distribution
    [p_group.muPhi,p_group.SigmaPhi,p_group.a_vPhi,p_group.b_vPhi] = ...
            MFX_VBupdate(...
            priors_group.muPhi,...
            iV_phi,...
            mphi,...
            Vphi,...
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
            mtheta,...
            Vtheta,...
            p_group.a_vTheta,...
            p_group.b_vTheta,...
            priors_group.a_vTheta,...
            priors_group.b_vTheta,...
            ind.theta_ffx,...
            ind.theta_in);

    % calculate free energy
    F = MFX_F(D_fit,p_group,priors_group,ind)

    %figures
    %     o_group.F = F;
    %     o_group.it = it;
    %     o_group.ind = ind;
    %     [o_group.options] = VBA_displayMFX(p_sub,o_sub,p_group,o_group);

    %figures
    %     figure(fig)
    %     hold on
    %     subplot(1,2,1);
    %     scatter(it,F(it));
    %     title('F')
    %     drawnow

    % save
    %     save(fullfile(S.path.hgf,['CORE_fittedparameters_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '_it' num2str(it) '.mat']),'D_fit','S','mphi','Vphi','mtheta','Vtheta','p_group','F');


    % stop rule
    %     if it>1
    %         dF = F(it) - F(it-1)
    % %         figure(fig)
    % %         scatter(it,abs(dF));
    % %         title('abs dF')
    %         
    %         if abs(dF) <= opt.TolFun || it >= opt.MaxIter
    %             stop = 1;
    %         end
    %     end
    %     
    %     it = it +1;

    % update priors in S
    r.c_obs.priormus=p_group.muPhi';
    r.c_obs.priorsas = transpose(p_group.b_vPhi./p_group.a_vPhi);
    S.c_obs = r.c_obs;
    r.c_prc.priormus=p_group.muTheta';
    r.c_prc.priorsas = transpose(p_group.b_vTheta./p_group.a_vTheta);
    S.c_prc = r.c_prc; 
    S.p_group = p_group;
    if ~isfield(S,'F')
        S.F=[];
    end
    S.F = [S.F F];

    % Update input files for next iteration
    in_dir = fullfile(pwd,'Data/input_files');
    for d = 1:nsub
        
        % input file index
        ii = nrand*(nr-1) + (d-1);
        
        % write variable
        in_fname = fullfile(in_dir,['input' num2str(ii) '.mat']);
        m = matfile(in_fname,'Writable',true);
        m.S = S;
    end

    % cleanup
    for d = 1:nsub
        
        % input file index
        ii = nrand*(nr-1) + (d-1);
        
        % delete
        fname = ['output' num2str(ii) '.mat'];
        disp(['copying ' fname])
        copyfile(fname,in_dir)
        delete(fname);
        delete(['input' num2str(ii) '.mat']);
    end
end
quit

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

function [in,out] = CORE_condor_monitor_outputs

nowtime = now;
in = dir('input*.mat');

fin=0;
while fin==0
    pause(10)
    out = dir('output*.mat');
    if length(out)==length(in)
        complete = [out(:).bytes]>0 & [out(:).datenum]>nowtime;
        if all(complete)
            fin=1;
        end
        disp(['number of outputs complete: ' num2str(sum(complete)) '/' num2str(length(complete))])
    end
end