function DCM_PEB_trim(GCMfile,DCMdir,run_num,fullmodel,loadorsave)
dbstop if error
disp('loading GCM...')
load(fullfile(DCMdir,GCMfile));
disp('finished loading GCM')
cd(DCMdir)
DCM = GCM{1,fullmodel};
if isempty(DCM.Bc)
    fields = {'A','C'};
else
    fields = {'B','C'};
end

% find out if a second level model is needed
if length(unique(Xb))>1
    seclev=1;
end

GCMex=GCM;
if excl_full
    GCMex(:,fullmodel) =[];
end
% first level model
%--------------------------------------------------------------------------
% Bayesian model averages, weighting each model by its marginal likelihood
% Applied over subjects using FFX Baysian parameter averaging
%try; bma  = spm_dcm_bma(GCMex);end;

% second level model
%--------------------------------------------------------------------------
M    = struct('X',Xb);

% extract results
[Ns Nm] = size(GCMex);
clear P Q R
Y=[];
Q=[];
F=[];
%for i = 1:Ns
%    
%    % data - over subjects
%    %----------------------------------------------------------------------
%    for nt = 1:length(GCM{i,1}.xY.y)
%        Y(:,i,nt) = GCM{i,1}.xY.y{nt}*DCM.M.U(:,1);
%    end
%    % Parameter averages (over models)
%    %----------------------------------------------------------------------
%    try; Q(:,i,1) = full(spm_vec(bma.SUB(i).Ep));end;
%    
%    % Free energies
%    %----------------------------------------------------------------------
%    for j = 1:Nm
%        F(i,j,1) = GCM{i,j}.F - GCM{i,1}.F;
%    end
%end
pst = GCM{1,1}.xY.pst;
clear GCMex bma

% Bayesian model reduction (avoiding local minima over models)
%==========================================================================
% spm_dcm_bmr operates on different DCMs of the same data (rows) to find
% the best model. It assumes the full model - whose free-parameters are
% the union (superset) of all free parameters in each model - has been
% inverted, and from that model creates reduced versions for the other
% models. If the other models have already been inverted, these are ignored
% and only the full model inversion is used to derive the others.
% GCM is a {Nsub x Nmodel} cell array of DCM filenames or model structures  
%         of Nsub subjects, where each model is reduced independently
disp('creating/loading RCM...')
rname=['RCM_fit' num2str(run_num) '.mat'];
if loadorsave && exist(fullfile(DCMdir,rname),'file')
    clear GCM
    load(fullfile(DCMdir,rname));
    try; load(fullfile(DCMdir,'data_and_estimates.mat'));end
else
    RCM   = spm_dcm_bmr(GCM,fields);
    clear GCM
    if iscell(RCM{1,1}); RCM = vertcat(RCM{:});end
    if excl_full
        RCM(:,fullmodel) =[];
    end
    % save
    try 
        save(rname,'RCM','BMC','xp','-v7.3');
    catch
        save(rname,'RCM','-v7.3');
    end
end

disp('finished creating/loading RCM')
% Bayesian model averages, weighting each model by its marginal likelihood pooled over subjects
%rma  = spm_dcm_bma(RCM);

if seclev
    % OPTION 1: Bayesian model reduction over the *joint* space of first and second level models
    % this empirical Bayesian approach is implemented so that one can compare different second level models without having to repeat the inversion of each subject's DCM. 
    % This can be thought of as a generalization of the standard summary statistic approach.

    % BMC - search over first and second level effects
    %--------------------------------------------------------------------------
    % This Bayesian model comparison should be contrasted with model
    % comparison at the second level. Here, we are interested in the best model
    % of first level parameters that show a second level effect. This is not
    % the same as trying to find the best model of second level effects (see next section for that). Model
    % comparison among second level parameters uses spm_dcm_peb_bmc.
    [BMC,PEB] = spm_dcm_bmc_peb(RCM,M,fields);

    % BMA - exhaustive search over second level parameters
    %--------------------------------------------------------------------------
    % In either case (spm_dcm_peb or spm_dcm_bmc_peb), one can now apply Bayesian model reduction to the posterior densities 
    % over the second level parameters (spm_dcm_peb_bmc.m) to find out where the key between-subject effects are expressed; 
    % in other words, identify the parameters that mediate group effects.
    % search over nested models
    PEB.gamma = 1/128;
    BMA       = spm_dcm_peb_bmc(PEB);
end

% posterior predictive density and LOO cross validation
%==========================================================================
%if length(unique(Xb))>1
%    spm_dcm_loo(RCM(:,1),Xb,fields);
%end
for i = 1:Ns
    % Parameter averages (over models)
    %----------------------------------------------------------------------%

    try; Q(:,i,2) = full(spm_vec(rma.SUB(i).Ep));end

    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j,2) = RCM{i,j}.F - RCM{i,1}.F;
    end

end
%[~,~,xp] = spm_dcm_bmc(RCM);
clear rma PEB

% hierarchical (empirical Bayes) model reduction:
% optimises the empirical priors over the parameters of a set of first level DCMs, using second level or
% between subject constraints specified in the design matrix X.
% See: https://en.wikibooks.org/wiki/SPM/Parametric_Empirical_Bayes_(PEB)
%==========================================================================
%[PEB,P]   = spm_dcm_peb(P,M,field)
% 'Field' refers to the parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
% 'All' will invoke all fields. This argument effectively allows 
%  one to specify the parameters that constitute random effects. 
% If P is an an {N [x M]} array, will return one PEB for all
% models/subjects
% If P is an an (N x M} array, will return a number of PEBs as a 
% {1 x M} cell array.
disp('creating/loading PCM...')
pname=['PCM_peb' num2str(run_num) '.mat'];
if loadorsave && exist(fullfile(DCMdir,pname),'file')
    clear RCM
    load(fullfile(DCMdir,pname));
    if isempty(F);try;load(fullfile(DCMdir,'data_and_estimates.mat'));end;end
else
    
    
    % OPTION 2: 
    % instead of empirical Bayesian model reduction considering every first level model (using spm_dcm_bmc_peb.m), 
    % instead we could just estimate the second level effects under the full model using (spm_dcm_peb.m).
    % Generates second level posteriors for each model that generally correspond to the group means and differences. 
    % These constitute empirical priors that shrink subject-wise estimates, thereby eliminating a degree of between subject variability. 
    [peb,PCM] = spm_dcm_peb(RCM,[],fields);
    clear RCM
    save(pname,'PCM','-v7.3');
end
% if 
%if strcmp(fields,'B')
    % In either case (spm_dcm_peb or spm_dcm_bmc_peb), one can now apply Bayesian model reduction to the posterior densities 
    % over the second level parameters (spm_dcm_peb_bmc.m) to find out where the key between-subject effects are expressed; 
    % in other words, identify the parameters that mediate group effects.
    % search over nested models
%    pma = spm_dcm_peb_bmc(peb);

    % Review results
    %spm_dcm_peb_review(pma,DCM);

    % then averages over models for classical 2nd level inference
    %pma  = spm_dcm_bma(PCM);
%end

for i = 1:Ns
    % Parameter averages (over models)
    %----------------------------------------------------------------------
    try;Q(:,i,3) = full(spm_vec(pma.SUB(i).Ep));end

    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j,3) = PCM{i,j}.F - PCM{i,1}.F;
    end

end

clear PCM
%clear PCM peb pma
%disp('finished creating/loading PCM')
save('data_and_estimates.mat','Y','Q','F','Xb','pst')

DCM_PEB_PLOT('data_and_estimates.mat')
