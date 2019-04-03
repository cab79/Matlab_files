function varargout=HGF_group_model_comparison(S,varargin)
% OUTPUT: varargout is {posterior,out,gposterior,gout}
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

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
                try
                    fname=[S.fname_pref '_percmodel' num2str(S.perc_models(pm)) '_respmodel' num2str(S.resp_models(rm)) S.fname_ext{om}];
                    ls = load(fullfile(fitted_path,fname));
                catch
                    fname=[S.fname_pref '_pm' num2str(S.perc_models(pm)) '_rm' num2str(S.resp_models(rm)) S.fname_ext{om}];
                    ls = load(fullfile(fitted_path,fname));
                end
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
%         tapas_fit_plotCorr_group(rs(mi,:));
%         title(['model ' num2str(mi) ', simple average'])
        try 
            rc_out=tapas_bayesian_parameter_average_CAB(1,rc(mi,:));
%             tapas_fit_plotCorr(rc_out);
%             title(['model ' num2str(mi) ', bayesian average'])
        end
    else
        rc_out = [];
    end
end
varargout={LMEgrp,rc_out};

% plot mean LME; separate by groups
if (rm>1 || pm>1) && om==1 
    
    % pm on x axis
%     pLME = reshape(sum(cat(2,LMEgrp{:}),2),rm,pm);
%     clims = [min(pLME(:)), max(pLME(:))];
%     figure
%     ax1=subplot('Position',[0.1 0.2 0.3 0.6])
%     imagesc(pLME,clims)
%     ax1.XAxisLocation = 'top';
%     xlabel('perceptual models')
%     ylabel('response models')
%     title(['Log-model evidence'])
%     ax2=subplot('Position',[0.42 0.2 0.1 0.6])
%     imagesc(mean(pLME,2),clims)
%     set(gca,'XTick',1,'XTickLabel','mean')
%     set(gca,'YTick',[])
% %     ax2.YAxisLocation = 'right';
%     ax2.XAxisLocation = 'top';
%     ax3=subplot('Position',[0.1 0.08 0.3 0.1])
%     imagesc(mean(pLME,1),clims)
%     set(gca,'YTick',1,'YTickLabel','mean')
%     set(gca,'XTick',[])
%     colorbar('Position',[0.54 0.08 0.03 0.72])
    
    pLME = reshape(sum(cat(2,LMEgrp{:}),2),rm,pm)';
    clims = [min(pLME(:)), max(pLME(:))];
    image_plot(pLME,S,clims,'Log-model evidence','mean',{})
    
    % plot mean LME; separate by groups
%     for g = 1:length(grpuni)
%         % mean over subjects
%         pLME = reshape(sum(LMEgrp{1,g},2),rm,pm);
%         figure
%         imagesc(pLME)
%         colorbar
%         xlabel('perc models')
%         ylabel('resp models')
%         title(['sum LME, group ' num2str(g)])
%     end
end

% compare models
% if length(LMEgrp)==1
%     [posterior,out] = VBA_groupBMC(LME);
%     varargout(1:2) = {posterior,out};
% elseif length(LMEgrp)>1
%     [gposterior,gout] = VBA_groupBMC_btwGroups_CAB(LMEgrp)
%     varargout(3:4) = {gposterior,gout};
% end

options.DisplayWin = 0;
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

% plot Efs - model frequencies
ef = pm_out{1}.Ef';
ef_pm = pm_out{1}.families.Ef';
ef_rm = rm_out{1}.families.Ef';
ef = reshape(ef,rm,pm)';
clims = [min(ef(:)), max(ef(:))];
image_plot(ef,S,clims,'Estimated model frequencies','family',{ef_pm,ef_rm})
% figure
% ax1=subplot('Position',[0.1 0.2 0.3 0.6])
% imagesc(ef,clims)
% ax1.XAxisLocation = 'top';
% xlabel('perceptual models')
% ylabel('response models')
% title(['Estimated model frequencies'])
% ax2=subplot('Position',[0.42 0.2 0.1 0.6])
% imagesc(ef_rm',clims)
% set(gca,'XTick',1,'XTickLabel','family')
% set(gca,'YTick',[])
% %     ax2.YAxisLocation = 'right';
% ax2.XAxisLocation = 'top';
% ax3=subplot('Position',[0.1 0.08 0.3 0.1])
% imagesc(ef_pm,clims)
% set(gca,'YTick',1,'YTickLabel','family')
% set(gca,'XTick',[])
% colorbar('Position',[0.54 0.08 0.03 0.72])

% plot EPs
ep = pm_out{1}.ep;
ep_pm = pm_out{1}.families.ep;
ep_rm = rm_out{1}.families.ep;
ep = reshape(ep,rm,pm)';
clims = [0 1];
image_plot(ep,S,clims,'Exceedence probability (EP)','family',{ep_pm,ep_rm})
% figure
% ax1=subplot('Position',[0.1 0.2 0.3 0.6])
% imagesc(ep,clims)
% ax1.XAxisLocation = 'top';
% xlabel('perceptual models')
% ylabel('response models')
% title(['Exceedence probability (EP)'])
% ax2=subplot('Position',[0.42 0.2 0.1 0.6])
% imagesc(ep_rm',clims)
% set(gca,'XTick',1,'XTickLabel','family')
% set(gca,'YTick',[])
% %     ax2.YAxisLocation = 'right';
% ax2.XAxisLocation = 'top';
% ax3=subplot('Position',[0.1 0.08 0.3 0.1])
% imagesc(ep_pm,clims)
% set(gca,'YTick',1,'YTickLabel','family')
% set(gca,'XTick',[])
% colorbar('Position',[0.54 0.08 0.03 0.72])

% plot PEPs
pep = pm_out{1}.pep;
pep_pm = pm_out{1}.families.pep;
pep_rm = rm_out{1}.families.pep;
pep = reshape(pep,rm,pm)';
clims = [0 1];
image_plot(pep,S,clims,'Protected exceedence probability (PEP)','family',{pep_pm,pep_rm})
% figure
% ax1=subplot('Position',[0.1 0.2 0.3 0.6])
% imagesc(pep,clims)
% ax1.XAxisLocation = 'top';
% xlabel('perceptual models')
% ylabel('response models')
% title(['Protected exceedence probability (PEP)'])
% ax2=subplot('Position',[0.42 0.2 0.1 0.6])
% imagesc(pep_rm',clims)
% set(gca,'XTick',1,'XTickLabel','family')
% set(gca,'YTick',[])
% %     ax2.YAxisLocation = 'right';
% ax2.XAxisLocation = 'top';
% ax3=subplot('Position',[0.1 0.08 0.3 0.1])
% imagesc(pep_pm,clims)
% set(gca,'YTick',1,'YTickLabel','family')
% set(gca,'XTick',[])
% colorbar('Position',[0.54 0.08 0.03 0.72])

% plot num of unique priors
% figure
% for mi = 1:size(priors,3)
%     for pi=1:2 % om priors
%         subplot(size(priors,3))
%         [~,~,ui]=length(unique(priors(mi,:,pi)));
%     end
% end
function image_plot(data,S,clims,title_text,XYlabel,family)
if isempty(family)
    family{1} = mean(data,2)';
    family{2} = mean(data,1);
end
figure
ax1=subplot('Position',[0.1 0.2 0.5 0.3]);
imagesc(data,clims)
set(gca,'XTick',S.resp_models,'YTick',S.perc_models,'FontSize',12)
ax1.XAxisLocation = 'top';
ylabel('perceptual models')
xlabel('response models')
title(title_text)
ax2=subplot('Position',[0.62 0.2 0.5/6 0.3])
imagesc(family{1}',clims)
set(gca,'XTick',1,'XTickLabel',XYlabel,'FontSize',12)
set(gca,'YTick',[])
%     ax2.YAxisLocation = 'right';
ax2.XAxisLocation = 'top';
ax3=subplot('Position',[0.1 0.08 0.5 0.1])
imagesc(family{2},clims)
set(gca,'YTick',1,'YTickLabel',XYlabel,'FontSize',12)
set(gca,'XTick',[])
colorbar('Position',[0.74 0.08 0.03 0.42])

colormap(flipud(cbrewer('seq', 'Reds', 100, 'pchip')))