function varargout=HGF_group_model_comparison(S)
% OUTPUT: varargout is {posterior,out,gposterior,gout}

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
pm_family{1,length(S.perc_models)} = [];
rm_family{1,length(S.resp_models)} = [];
for pm = 1:length(S.perc_models)
%     load(fullfile(fitted_path,['CORE_fittedparameters_percmodel' num2str(S.perc_models(m)) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']));

    % testing response models
    for rm = 1:length(S.resp_models)
        mi=mi+1;   
        
        % factor matrix and families
        fm(mi,1) = pm;
        fm(mi,2) = rm;
        pm_family{pm} = [pm_family{pm} mi];
        rm_family{rm} = [rm_family{rm} mi];
        
        fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext];
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
        end
    end
end

% remove bad subjects
mr = unique(mr);
LME(:,mr) = [];
rs(:,mr) = [];
rc(:,mr) = [];
grplist(mr) = [];
    
for mi = 1:size(LME,1)
    
    % separate into groups
    for g = 1:length(grpuni)
        LMEgrp{1,g}(mi,:) = LME(mi,strcmp(grplist,grpuni{g}));
    end

    % diagnostics: parameter correlations
    tapas_fit_plotCorr_group(rs(mi,:));
    title(['model ' num2str(mi) ', simple average'])
    try 
        rc_out=tapas_bayesian_parameter_average_CAB(1,rc(mi,:));
        tapas_fit_plotCorr(rc_out);
        title(['model ' num2str(mi) ', bayesian average'])
    end
end

% plot mean LME; separate by groups
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

% compare models
% if length(LMEgrp)==1
%     [posterior,out] = VBA_groupBMC(LME);
%     varargout(1:2) = {posterior,out};
% elseif length(LMEgrp)>1
%     [gposterior,gout] = VBA_groupBMC_btwGroups_CAB(LMEgrp)
%     varargout(3:4) = {gposterior,gout};
% end

% compare perc model families
options.families = pm_family;
if length(LMEgrp)==1
    [~,~,pm_out] = VBA_groupBMC(LME,options);
elseif length(LMEgrp)>1
    [~,~,pm_out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options)
end

% compare resp model families
options.families = rm_family;
if length(LMEgrp)==1
    [~,~,rm_out] = VBA_groupBMC(LME,options);
elseif length(LMEgrp)>1
    [~,~,rm_out] = VBA_groupBMC_btwGroups_CAB(LMEgrp,options)
end
varargout={pm_out,rm_out};

