clear all
close all
models = {
    %'ExpStudy_HGF_results_Sm5';
    'ExpStudy_HGF_results_Sm10';
    'ExpStudy_HGF_results_Sm20';
    'ExpStudy_HGF_results_Sm40';
    'ExpStudy_HGF_results_LLS1';
    %'ExpStudy_HGF_results_Sm5';
    %'ExpStudy_HGF_results_Sm6';
    %'ExpStudy_HGF_results_Sm7';
    %'ExpStudy_HGF_results_Sm8';
    %'ExpStudy_HGF_results_smNoc_RW-1';
    %'ExpStudy_HGF_results_Sm5_RW-1';
    %'ExpStudy_HGF_results_RW-1';
    %'ExpStudy_HGF_results_v4';
    };

for m = 1:length(models)
    fname = [models{m} '.xlsx'];
    T = readtable([fname]);

    g1(m,:) = T.LME_av(1:16)';
    g2(m,:) = T.LME_av(17:31)';
    g3(m,:) = T.LME_av(32:47)';
    g_all(m,:) = T.LME_av';
end

[h, p] = VBA_groupBMC_btwGroups({g1, g2, g3})
%[g1_posterior,g1_out] = VBA_groupBMC(g1);
%[g2_posterior,g2_out] = VBA_groupBMC(g2);
%[g_all_posterior,g_all_out] = VBA_groupBMC(g_all);

