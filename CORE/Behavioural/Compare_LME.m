clear all
close all

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
    
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_mmtrials\results';
    
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd05_0-150\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_mm_thresh_b0_sd05_0-150\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm_thresh_b0_sd05_0-150\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_mm_thresh_b0_sd05_0-150\results';
    
    };

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
pdatfile = 'C:\Data\CORE\Participant_data.xlsx';
% names of headers in the above xls file:
    grphead = 'Group';
    inchead = 'Include';
% which codes to analyse in 'Include' columns in participant data file?
include_codes = [1];
[~,~,pdata] = xlsread(pdatfile);
grp_col = find(strcmp(pdata(1,:),grphead));
inc_col = find(strcmp(pdata(1,:),inchead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));
grp = pdata(2:end,grp_col);
[grpuni,~,grpnum] = unique(grp);
grpnum = grpnum(inc_idx);

%% Compare LME
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


if 0
    %% Compare parameters
    paraname = 'al0'; % first part of name common to all names
    % compare rows (meaned within rows)
    paranum = [1;2]; % 
               %3;4]; % 
    for m = 1:length(models)
        fname = [models{m} '.xlsx'];
        T = readtable([fname]);
        hdr = T.Properties.VariableNames;
        temp = strfind(hdr,paraname);
        hi = find(not(cellfun('isempty', temp)));
        xy=[];
        for pn = 1:size(paranum,1)
            pi=paranum(pn,:);
            pv=[];
            for pii = 1:length(pi)
                pv(:,pii)=T.(hdr{hi(paranum(pn,pii))});
            end
            xy(:,pn)=reshape(pv,1,[]);
        end
        pval(m) = friedman(xy,pii,'off');
        %figure;errorbar(1:size(xy,2),mean(xy),std(xy))
        cind = repmat(1:size(xy,2),size(xy,1),1);
        figure;line(1:size(xy,2),xy);
        title(models{m})
    end
end