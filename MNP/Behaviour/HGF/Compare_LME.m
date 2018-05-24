clear all
close all

compareLME=1;
compareParam=0;

addpath(genpath('C:\Data\Matlab\VBA-toolbox-master'));

models = {
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\Alpha7\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\Alpha4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\Alpha2a\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\Alpha2b\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\Alpha1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-2alpha\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_varPrior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_varP_noRB\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_varP_noRB_1alpha\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_noRB_1alpha\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-4priors\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-4alphas\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-4alpri\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF2\results'; % parameter estimability is poor 
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_soft_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_soft\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_soft\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_soft\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RT_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RT\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_RT\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_RT\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RTsoft_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RTsoft\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_RTsoft\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_RTsoft\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_soft_a4_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_soft_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_soft_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_soft_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RT_a4_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RT_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_RT_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_RT_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RTsoft_a4_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RTsoft_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_RTsoft_a4\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_RTsoft_a4\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_soft_a2_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_soft_a2\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_soft_a2\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_soft_a2\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RT_a2_noprior\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_RT_a2\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_RT_a2\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_RT_a2\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_sd1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_sd1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_sd1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_sd1\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd1\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd05\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd05\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd05\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd05\results';
    
   % 'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd0\results';
   % 'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd0\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd0\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd0\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd15\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd15\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd15\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd15\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_mmtrials\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd05_0-150\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd05_0-150\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd05_0-150\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd05_0-150\results';
    
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd05_0-150_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd05_0-150_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd05_0-150_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd05_0-150_mmtrials\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mmERPrun1_thresh_b0_sd05_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mmERPrun1_thresh_b0_sd05_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mmERPrun1_thresh_b0_sd05_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mmERPrun1_thresh_b0_sd05_mmtrials\results';
    
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mmERPrun1_thresh_b0_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mmERPrun1_thresh_b0_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mmERPrun1_thresh_b0_mmtrials\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mmERPrun1_thresh_b0_mmtrials\results';
    };

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
pdatfile = 'C:\Data\CORE\Participant_data.xlsx';
% names of headers in the above xls file:
    grphead = 'Group';
    inchead = 'Include';
    dathead = {'Acc_aff'}; % can be all of them to be analysed, or just the first one (see next line)
    alldat_head = 1; % analysis all other dat columns too
    
% which codes to analyse in 'Include' columns in participant data file?
include_codes = [1];
[~,~,pdata] = xlsread(pdatfile);
grp_col = find(strcmp(pdata(1,:),grphead));
inc_col = find(strcmp(pdata(1,:),inchead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));
grp = pdata(2:end,grp_col);
nn = cell2mat(cellfun(@(x) ischar(x),grp,'uniformoutput',0))
[grpuni,~,grpnum] = unique(grp(nn));
grpnum = grpnum(inc_idx);

% extract data for correlation with model parameters
dat_col = find(ismember(pdata(1,:),dathead));
if alldat_head
    dat_col = dat_col:size(pdata,2);
end
dat = pdata(2:end,dat_col);
dat = cell2mat(dat(grpnum>0,:));
dat_head = pdata(1,dat_col);

%% Compare LME
if compareLME
    for m = 1:length(models)
        fname = [models{m} '.xlsx'];
        T = readtable([fname]);

        for g = 1:length(grpuni)
            gdat{1,g}(m,:) = T.LME_av(grpnum==g)';
        end
        g_all(m,:) = T.LME_av(grpnum>0)';
    end
    [h, p] = VBA_groupBMC_btwGroups(gdat)
    %[g_all_posterior,g_all_out] = VBA_groupBMC(g_all);
end


if compareParam
    %% Compare parameters
    paraname = {'al0','al1','rb','om','dau','da','ud','wt','psi','epsi','mu','sa','muhat','sahat','be0','be1','be2','be3','be4','be5','ze','be'}; % first part of name common to all names
    % compare rows (meaned within rows)
    paranum = [1,2]; % 
               %3;4]; % 
    paras={};
    paracorr=[]
    for m = 1:length(models)
        fname = [models{m} '.xlsx'];
        T = readtable([fname]);
        hdr = T.Properties.VariableNames;
        for pm = 1:length(paraname)
            temp = strfind(hdr,paraname{pm});
            hi = find(not(cellfun('isempty', temp)));
            xy=[];
            if ~isempty(hi)
                for pn = 1:size(paranum,1)
                    pi=paranum(pn,:);
                    pv=[];
                    for pii = 1:min(length(pi),length(hi))
                        temp = T.(hdr{hi(paranum(pn,pii))});
                        pv(:,pii)=temp(grpnum>0);
                    end
                    %xy(:,pn)=reshape(pv,1,[]);
                    xy(:,pn)=mean(pv,2);
                end
                paras{m,pm} = xy;
            else
                paras{m,pm} = nan(length(grpnum>0),1);
            end
                
            if size(paranum,1)>1
                pval(m) = friedman(xy,pii,'off');
                %figure;errorbar(1:size(xy,2),mean(xy),std(xy))
                cind = repmat(1:size(xy,2),size(xy,1),1);
                figure;line(1:size(xy,2),xy);
                %title(models{m})
                ax = gca;
                ax.XLim = [0.9 2.1];
                ax.XTick = [1 2];
                set(gca,'XTickLabel',{'DC1','DC3'})
                set(gca, 'FontSize', 20)
            end
        end
    end
    
    % calculate cross-correlations over model combinations
    %nc = combnk(1:m,2);
    %for pn = 1:size(paranum,1)
    %    for n = 1:length(nc)
    %        [~,paracorr(n,pn)] = corr(paras{nc(n,1)}(:,pn),paras{nc(n,2)}(:,pn),'type','Spearman');
    %    end
    %end
    
    % correlate with behavioural variables
    para_results = struct;
    para_results.dat_headers = dat_head;
    para_results.models = models;
    for pm = 1:length(paraname)
        for d = 1:size(dat,2)
            for m = [4 8]%1:length(models)
                x=paras{m,pm};%(:,pn);
                y=dat(:,d);
                rm=isnan(x)+isnan(y);
                x=x(rm==0);
                y=y(rm==0);
                if ~isempty(x) && ~isempty(y)
                    [~,paracorr(m,d)] = corr(x,y,'type','Spearman');
                    if paracorr(m,d)<0.01
                        figure
                        scatter(tiedrank(x),tiedrank(y)); title([paraname{pm} ', model ' num2str(m) ', ' dat_head{d}]);
                        %scatter(x,y); title([paraname{pm} ', model ' num2str(m) ', ' dat_head{d}]);
                    else
                        paracorr(m,d) = NaN;
                    end
                else
                    paracorr(m,d) = NaN;
                end
            end
        end
        para_results.(paraname{pm}) = paracorr;
    end
    
end