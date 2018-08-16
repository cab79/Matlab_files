function varargout=HGF_group_model_comparison(S)
% OUTPUT: varargout is {posterior,out,gposterior,gout}

fitted_path = fullfile(S.path.hgf,'fitted');
cd(fitted_path)

% group info
grplist = S.designmat(2:end,strcmp(S.designmat(1,:),'groups'));
grpuni = unique(grplist,'stable');

% for each model
rm=[];
rc=cell(1,1);

% testing perceptual models
% for m = 1:length(S.perc_models)
%     load(fullfile(fitted_path,['CORE_fittedparameters_percmodel' num2str(S.perc_models(m)) '_respmodel' num2str(S.resp_model) '_fractrain' num2str(S.frac_train) '_' S.sname '.mat']));

% testing response models
for m = 1:length(S.resp_models)
    ls = load(fullfile(fitted_path,[S.fname_pref S.resp_models{m} S.fname_ext]));
    D_fit=ls.D_fit;
    % mean evidence for each fitted model over repetitions
    for d = 1:length(D_fit)
        try
            LME(m,d) = D_fit(d).HGF.fit.optim.LME; % rows: fitted models; columns: simulated models
        catch
            rm=[rm d];
            try
                disp(['Missing LME from subject d = ' num2str(d) ', model ' num2str(S.perc_models(m))])
            end
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
    
%for m = 1:length(S.perc_models)
for m = 1:length(S.resp_models)
    
    % separate into groups
    for g = 1:length(grpuni)
        LMEgrp{1,g}(m,:) = LME(m,strcmp(grplist,grpuni{g}));
    end
    
    % diagnostics: parameter correlations
    tapas_fit_plotCorr_group(rs(m,:));
    rc_out=tapas_bayesian_parameter_average_CAB(1,rc(m,:));
    tapas_fit_plotCorr(rc_out);
end

% compare models
if length(LMEgrp)==1
    [posterior,out] = VBA_groupBMC(LME);
    varargout(1:2) = {posterior,out};
elseif length(LMEgrp)>1
    [gposterior,gout] = VBA_groupBMC_btwGroups(LMEgrp)
    varargout(3:4) = {gposterior,gout};
end
