function [posterior,out,gposterior,gout]=HGF_group_model_comparison(S)

fitted_path = fullfile(S.path.hgf,'fitted');
cd(fitted_path)

% group info
grplist = S.designmat(2:end,strcmp(S.designmat(1,:),'groups'));
grpuni = unique(grplist,'stable');

% for each model
rm=[];
rc=cell(1,1);
for m = 1:length(S.perc_models)
    load(fullfile(fitted_path,['CORE_fittedparameters_percmodel' num2str(S.perc_models(m)) '_respmodel' num2str(S.resp_model) '_' S.sname '.mat']));
    
    
    % mean evidence for each fitted model over repetitions
    for d = 1:length(D_fit)
        try
            LME(m,d) = D_fit(d).HGF.fit.optim.LME; % rows: fitted models; columns: simulated models
        catch
            rm=[rm d];
            disp(['Missing LME from subject d = ' num2str(d) ', model ' num2str(S.perc_models(m))])
        end
        rs(m,d)=D_fit(d).HGF.fit;
        rc{m,d}=D_fit(d).HGF.fit;
    end
end

% remove bad subjects
rm = unique(rm);
LME(:,rm) = [];
rs(:,rm) = [];
rc(:,rm) = [];
grplist(rm) = [];
    
for m = 1:length(S.perc_models)
    
    % separate into groups
    for g = 1:length(grpuni)
        LMEgrp{1,g}(m,:) = LME(m,strcmp(grplist,grpuni{g}));;
    end
    
    % diagnostics: parameter correlations
    tapas_fit_plotCorr_group(rs(m,:));
    rc_out=tapas_bayesian_parameter_average_CAB(1,rc(m,:));
    tapas_fit_plotCorr(rc_out);
end

% compare models
%[posterior,out] = VBA_groupBMC(LME);
[gposterior,gout] = VBA_groupBMC_btwGroups(LMEgrp)

