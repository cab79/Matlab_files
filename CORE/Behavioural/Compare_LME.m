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
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF2-alt\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev-soft\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-soft\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-soft\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev-RT\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RT\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RT-newprior\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-RT\results'; % parameter estimability is poor
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev-RTsoft\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_RTsoft_new\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-RTsoft\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft-b0\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft-b1\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft-b5\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft-b0s0\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft-b1s0\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF-RTsoft-b5s0\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-RTsoft-b0s0\results'; 
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-RTsoft-b1s0\results';
    %'C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT-RTsoft-b5s0\results'; 
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mm\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mmtrials\results';
    'C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mmtrials_mmpos\results';
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
grp = pdata(1+inc_idx,grp_col);
[grpuni,~,grpnum] = unique(grp);

for m = 1:length(models)
    fname = [models{m} '.xlsx'];
    T = readtable([fname]);

    for g = 1:length(grpuni)
        gdat{1,g}(m,:) = T.LME_av(grpnum==g)';
    end
    g_all(m,:) = T.LME_av';
end

[h, p] = VBA_groupBMC_btwGroups(gdat)
%[g_all_posterior,g_all_out] = VBA_groupBMC(g_all);

