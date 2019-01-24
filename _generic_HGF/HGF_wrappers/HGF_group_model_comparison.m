function varargout=HGF_group_model_comparison(S,varargin)
% OUTPUT: varargout is {posterior,out,gposterior,gout}

if nargin>1 && ~isempty(varargin{1})
    LME = varargin{1};
    LME_input = 1;
else
    LME_input = 0; % load from HGF structures
end

if ~isfield(S,'family_on')
    S.family_on=0;
end

if ~iscell(S.fname_ext)
    S.fname_ext = {S.fname_ext};
end

fitted_path = fullfile(S.path.hgf,'fitted');
cd(fitted_path)

% group info
grplist = S.designmat(2:end,strcmp(S.designmat(1,:),'groups'));
grpuni = unique(grplist,'stable');

% for each model
mr=[];
rc=cell(1,1);

% testing perceptual models
mi=0;
if S.family_on
    pm_family{1,length(S.perc_models)} = [];
    rm_family{1,length(S.resp_models)} = [];
end

for pm = 1:length(S.perc_models)
%     load(fullfile(fitted_path,['CORE_fittedparameters_percmodel' num2str(S.perc_models(m)) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']));

    % testing response models
    for rm = 1:length(S.resp_models)
        
        % same model, different runs
        for om = 1:length(S.fname_ext)
            mi=mi+1;   

            if S.family_on
                % factor matrix and families
                fm(mi,1) = pm;
                fm(mi,2) = rm;
                pm_family{pm} = [pm_family{pm} mi];
                rm_family{rm} = [rm_family{rm} mi];
            end

            if ~LME_input
                fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext{om}];
                ls = load(fullfile(fitted_path,fname));
                D_fit=ls.D_fit;
                % mean evidence for each fitted model over repetitions
                for d = 1:length(D_fit)
                    try
                        LME(mi,d) = D_fit(d).HGF.fit.optim.LME; % rows: fitted models; columns: simulated models
                    catch
                        mr=[mr d];
                        try
                            disp(['Missing LME from file: ' fname])
                        end
                    end
                    rs(mi,d)=D_fit(d).HGF.fit;
                    rc{mi,d}=D_fit(d).HGF.fit;

                    % find number of subjects with unique priors (due to inabilty to optimise at the commonly use prior)
                    priors(mi,d,:) = D_fit(d).HGF.fit.c_prc.PL.ommu(2:3);
                end
            end
        end
        
    end
end

if ~LME_input
    % remove bad subjects
    mr = unique(mr);
    LME(:,mr) = [];
    rs(:,mr) = [];
    rc(:,mr) = [];
    grplist(mr) = [];
    priors(:,mr,:) = [];
end
    
for mi = 1:size(LME,1)
    
    % separate into groups
    for g = 1:length(grpuni)
        LMEgrp{1,g}(mi,:) = LME(mi,strcmp(grplist,grpuni{g}));
    end

    if ~LME_input
        % diagnostics: parameter correlations
        tapas_fit_plotCorr_group(rs(mi,:));
        title(['model ' num2str(mi) ', simple average'])
        try 
            rc_out=tapas_bayesian_parameter_average_CAB(1,rc(mi,:));
            tapas_fit_plotCorr(rc_out);
            title(['model ' num2str(mi) ', bayesian average'])
        end
    else
        rc_out = [];
    end
end
varargout={LMEgrp,rc_out};

% plot mean LME; separate by groups
if (rm>1 || pm>1) && om==1 
    for g = 1:length(grpuni)
        % mean over subjects
        pLME = reshape(mean(LMEgrp{1,g},2),rm,pm);
        figure
        imagesc(pLME)
        colorbar
        xlabel('perc models')
        ylabel('resp models')
        title(['mean LME, group ' num2str(g)])
    end
end

% compare models
% if length(LMEgrp)==1
%     [posterior,out] = VBA_groupBMC(LME);
%     varargout(1:2) = {posterior,out};
% elseif length(LMEgrp)>1
%     [gposterior,gout] = VBA_groupBMC_btwGroups_CAB(LMEgrp)
%     varargout(3:4) = {gposterior,gout};
% end

if S.family_on
    % compare perc model families
    options.families = pm_family;
    if length(LMEgrp)==1
        [~,~,pm_out] = VBA_groupBMC_cab(LME,options,S.pep_flag);
    elseif length(LMEgrp)>1
        [~,~,pm_out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options,S.pep_flag);
    end

    % compare resp model families
    options.families = rm_family;
    if length(LMEgrp)==1
        [~,~,rm_out] = VBA_groupBMC_cab(LME,options,S.pep_flag);
    elseif length(LMEgrp)>1
        [~,~,rm_out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options,S.pep_flag);
    end
    varargout=[varargout {pm_out,rm_out}];
else
    options = {};
    % compare models
    if length(LMEgrp)==1
        [~,~,out] = VBA_groupBMC_cab(LME,options,S.pep_flag);
    elseif length(LMEgrp)>1
        [~,~,out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options,S.pep_flag);
    end
    varargout=[varargout {out}];
end

% plot num of unique priors
% figure
% for mi = 1:size(priors,3)
%     for pi=1:2 % om priors
%         subplot(size(priors,3))
%         [~,~,ui]=length(unique(priors(mi,:,pi)));
%     end
% end
