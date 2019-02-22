%% Analysis: Perceptual model fitting
% This function approaches model inversion from an empirical Bayes
% perspective (based on VBA_MFX function in VBA toolbox), 
% whereby within-subject priors are iteratively refined and
% matched to the inferred parent population distribution.
%  Note: all subjects must use the same model
% See: http://mbb-team.github.io/VBA-toolbox/wiki/VBA-MFX/

function CORE_condor_fit_step2

dbstop if error

pth = pwd;

% add toolbox paths
addpath(genpath(fullfile(pth, 'dependencies')))
addpath(genpath(fullfile(pth, 'Data')))

% get input and output file info
CORE_condor_monitor_outputs;

nsub = importdata('nsub.txt');
nrand = importdata('nrand.txt');
it = importdata('it.txt');
ngrp = length(nsub);

% nuisance covariate
datfile = fullfile(pth, 'Data', 'Participant_data.xlsx'); % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest
pdata = readtable(datfile);
Cov = {
   pdata.Age(find(pdata.Include));
    };

for nr = 1:nrand
    D_fit=struct;
    GS=struct;
    rstart = sum(nsub)*(nr-1); % starting index of randomisation
    
    
    for g = 1:ngrp
        gstart = (g-1)*nsub(max(g,2)-1); % starting index of group
        gend = gstart + nsub(g)-1; % ending index of group
        
        for d = 1:nsub(g)
            ns = gstart+d;

            % input file index
            ii = rstart + (ns-1);

            % create D_fit
            in_fname = ['input' num2str(ii) '.mat'];
            out_fname = ['output' num2str(ii) '.mat'];
            load(in_fname,'u','y');
            fit=load(out_fname,'out');
            D_fit(ns).HGF.u=u;
            D_fit(ns).HGF.y=y;
            D_fit(ns).HGF.fit=fit.out;
        end

        % get group prior stats
        ii = rstart + gstart; % index of first file in each group/rand
        load(['input' num2str(ii) '.mat'],'S');
        priors_group(g) = S.priors_group;
        p_group(g) = S.p_group;
        ind{g} = S.ind;
        iV_phi{g} = VBA_inv(priors_group(g).SigmaPhi);
        temp=priors_group(g).SigmaTheta;
        temp(isnan(temp)) = 0;
        iV_theta{g} = VBA_inv(temp); % uses temp because doesn't like nans

        % store sufficient statistics
        r(g).c_prc = eval([S.prc_config '(S,0)']); % runs the prc config function
        r(g).c_obs = eval([S.obs_config '(r(g),S)']); % runs the obs config function
        optim_prc_ind{g} = 1:length(ind{g}.theta_in);
        optim_obs_ind{g} = length(ind{g}.theta_in)+1:length(ind{g}.theta_in)+length(ind{g}.phi_in);
        for d = 1:nsub(g)
            ns = gstart+d;
            mphi{g}(:,d) = D_fit(ns).HGF.fit.p_obs.p;
            Vphi{g}{d} = nan(length(mphi{g}(:,d))); 
            Vphi{g}{d}(ind{g}.phi_in,ind{g}.phi_in) = D_fit(ns).HGF.fit.optim.Sigma(optim_obs_ind{g},optim_obs_ind{g});
            % for mtheta, need to log values outputted from HGF
            [pvec, ~] = GBM_transp_log(r(g),D_fit(ns).HGF.fit.p_prc.p);
            mtheta{g}(:,d) = pvec;
            % others are already ok as variances
            Vtheta{g}{d} = nan(length(mtheta{g}(:,d))); 
            Vtheta{g}{d}(ind{g}.theta_in,ind{g}.theta_in) = D_fit(ns).HGF.fit.optim.Sigma(optim_prc_ind{g},optim_prc_ind{g});
        end
        if g==1
            GS=S;
        else
            GS(g)=S;
        end
    end
        
    % correct for Cov
    % NEED TO ENSURE COV IS IN SAME ORDER AS SUBJECTS
    mphi_all = cat(2,mphi{:});
    mtheta_all = cat(2,mtheta{:});
    if ~isempty(Cov)
        covX=[ones(size(mphi_all,2),1) cell2mat(Cov)];
        for i = 1:size(mphi_all,1)
            if std(mphi_all(i,:))>0 && all(~isnan(mphi_all(i,:)))
                meanval = mean(mphi_all(i,:));
                [~,~,resid,~,~] = regress(mphi_all(i,:)',covX);
                mphi_all(i,:) = meanval + resid;
            end
        end
        for i = 1:size(mtheta_all,1)
            if std(mtheta_all(i,:))>0 && all(~isnan(mtheta_all(i,:)))
                meanval = mean(mtheta_all(i,:));
                [~,~,resid,~,~] = regress(mtheta_all(i,:)',covX);
                mtheta_all(i,:) = meanval + resid;
            end
        end
    end

    for g = 1:ngrp
%         S=GS(g);
        gstart = (g-1)*nsub(max(g,2)-1); % starting index of group
        gend = gstart + nsub(g)-1; % ending index of group
        
        mphi{g}=mphi_all(:,gstart+1:gend+1);
        mtheta{g}=mtheta_all(:,gstart+1:gend+1);
        
        % update moments of the parent population distribution
        [p_group(g).muPhi,p_group(g).SigmaPhi,p_group(g).a_vPhi,p_group(g).b_vPhi] = ...
                MFX_VBupdate(...
                priors_group(g).muPhi,...
                iV_phi{g},...
                mphi{g},...
                Vphi{g},...
                p_group(g).a_vPhi,...
                p_group(g).b_vPhi,...
                priors_group(g).a_vPhi,...
                priors_group(g).b_vPhi,...
                ind{g}.phi_ffx,...
                ind{g}.phi_in);
        [p_group(g).muTheta,p_group(g).SigmaTheta,p_group(g).a_vTheta,p_group(g).b_vTheta] = ...
                MFX_VBupdate(...
                priors_group(g).muTheta,...
                iV_theta{g},...
                mtheta{g},...
                Vtheta{g},...
                p_group(g).a_vTheta,...
                p_group(g).b_vTheta,...
                priors_group(g).a_vTheta,...
                priors_group(g).b_vTheta,...
                ind{g}.theta_ffx,...
                ind{g}.theta_in);

        % calculate free energy
        F = MFX_F(D_fit(gstart+1:gend+1),p_group(g),priors_group(g),ind{g})

        %figures
        %     o_group.F = F;
        %     o_group.it = it;
        %     o_group.ind{g} = ind{g};
        %     [o_group.options] = VBA_displayMFX(p_sub,o_sub,p_group(g),o_group);

        %figures
        %     figure(fig)
        %     hold on
        %     subplot(1,2,1);
        %     scatter(it,F(it));
        %     title('F')
        %     drawnow

        % save
        %     save(fullfile(S.path.hgf,['CORE_fittedparameters_percmodel' num2str(S.perc_model) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '_it' num2str(it) '.mat']),'D_fit','S','mphi{g}','Vphi{g}','mtheta{g}','Vtheta{g}','p_group(g)','F');


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
        r(g).c_obs.priormus=p_group(g).muPhi';
        r(g).c_obs.priorsas = transpose(p_group(g).b_vPhi./p_group(g).a_vPhi);
        GS(g).c_obs = r(g).c_obs;
        r(g).c_prc.priormus=p_group(g).muTheta';
        r(g).c_prc.priorsas = transpose(p_group(g).b_vTheta./p_group(g).a_vTheta);
        GS(g).c_prc = r(g).c_prc; 
        GS(g).p_group = p_group(g);
        if ~isfield(GS(g),'F')
            GS(g).F=[];
        end
        GS(g).F = [GS(g).F F];

        % Update input files for next iteration
        in_dir = fullfile(pwd,'Data/input_files');
        for d = 1:nsub(g)
            ns = gstart+d;

            % input file ind{g}ex
            ii = rstart + (ns-1);

            % write variable
            in_fname = fullfile(in_dir,['input' num2str(ii) '.mat']);
            m = matfile(in_fname,'Writable',true);
            m.S = GS(g);
        end
        
        %if g==1
        %    GS=S;
        %else
%             GS(g)=S;
        %end

        % cleanup
        for d = 1:nsub(g)
            ns = gstart+d;

            % input file ind{g}ex
            ii = rstart + (ns-1);

            % delete
            fname = ['output' num2str(ii) '.mat'];
            disp(['copying ' fname])
            copyfile(fname,in_dir)
            delete(fname);
            delete(['input' num2str(ii) '.mat']);
        end
    end
    
    % save D_fit,S
    save_dir = fullfile(pwd,'Data/fitted');
    if ~exist(save_dir)
        mkdir(save_dir)
    end
    save(fullfile(save_dir,['D_fit_r' num2str(nr) '_it' num2str(it)]),'D_fit','GS')
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
    
    % outputs
    out = dir('output*.mat');
    if length(out)==length(in)
        complete = [out(:).bytes]>0 & [out(:).datenum]>nowtime;
        if all(complete)
            fin=1;
        end
        disp(['number of outputs complete: ' num2str(sum(complete)) '/' num2str(length(complete))])
    end
    
    
    % errors
    out = dir('*.err');
    if ~isempty(out)
        fsize = [out(:).bytes];
        err_ind = find(fsize>0);
        for n = 1:length(err_ind)
            dat=importdata(out(err_ind(n)).name);
            dat(~cellfun(@isempty,dat))
        end
    end
end


function which_stragglers
% run this code manually
stragglers = [153]+1;
nmods = [1 6];
nsub = 44;
strag=struct;
i=0;
for pm = 1:nmods(1)
    for rm = 1:nmods(2)
        for d = 1:nsub
            % input file index
            ii = (pm-1)*nmods(2)*nsub + (rm-1)*nsub + d;
            if any(stragglers==ii)
                i=i+1;
                strag(i).pm = pm;
                strag(i).rm = rm;
                strag(i).d = d;
            end
        end
    end
end