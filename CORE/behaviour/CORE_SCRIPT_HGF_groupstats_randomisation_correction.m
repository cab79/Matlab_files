% extract, tabulate and save parameters and summary stats of trajectories

close all
clear all
dbstop if error
restoredefaultpath

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

% file info
S.path.hgf = 'C:\Data\CORE\behaviour\hgf'; 
S.path.hgf_rand = 'C:\Data\CORE\behaviour\hgf\fitted\randomised'; 

fnames = {
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it1.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it2.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it3.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it4.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it5.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it6.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it7.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it8.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it9.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it10.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it11.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it12.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it13.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it14.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it15.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it16.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it17.mat'
'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it18.mat'
};
rand_fnames = {
'D_fit_r*_it1.mat'
'D_fit_r*_it2.mat'
'D_fit_r*_it3.mat'
'D_fit_r*_it4.mat'
'D_fit_r*_it5.mat'
'D_fit_r*_it6.mat'
'D_fit_r*_it7.mat'
'D_fit_r*_it8.mat'
'D_fit_r*_it9.mat'
'D_fit_r*_it10.mat'
'D_fit_r*_it11.mat'
'D_fit_r*_it12.mat'
'D_fit_r*_it13.mat'
'D_fit_r*_it14.mat'
'D_fit_r*_it15.mat'
'D_fit_r*_it16.mat'
'D_fit_r*_it17.mat'
'D_fit_r*_it18.mat'
};

    
% which trial-wise stats to output for traj? Averages over trials
% without weighting by condition
S.summary_stats = {};%{'trial_mean','trial_std','trial_absmean','trial_absstd'}; % for all traj

% which traj to outputs to mean over conditions? First averages for
% each condition, before creating further averages between conditions
%S.condmean = {'PL_muhat_1','PL_dau','PL_sahat_1','PL_muhat_2','PL_sahat_2','PL_muhat_3','PL_sahat_3'};
%S.condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3','PL_psi_1','PL_psi_2','PL_psi_3','PL_sa_1','PL_sa_2','PL_sa_3','PL_sahat_1','PL_sahat_2','PL_sahat_3','PL_mu_1','PL_mu_2','PL_mu_3','PL_muhat_1','PL_muhat_2','PL_muhat_3'};
S.condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3'};

%% marginal means of interest
S.cond.condmean = [1:24]; %0.5 prob, 3 digit change, right hand

ranksum_all_p = struct;
ranksum_all_z = struct;
medians_all = struct;
for f=1:length(fnames)

    temp=load(fullfile(S.path.hgf,'fitted',fnames{f}));
    LS = temp.S;
    D_fit = temp.D_fit;
    
    frand = dir(fullfile(S.path.hgf_rand,rand_fnames{f}));
    D_rand=struct;
    for fr = 1:length(frand)
        disp(['loading HGF, randomisation ' num2str(fr)])
        temp=load(fullfile(S.path.hgf_rand,frand(fr).name));
        D_rand(fr).D_fit = temp.D_fit;
        D_rand(fr).S = temp.GS;
    end
    
    %% extract data
    S.designmat = LS.designmat;
    [out.T,out.traj,out.param,out.rt] = CORE_extract_HGF_results(D_fit,S);
    
    %% extract rand data
    for fr = 1:length(frand)
        disp(['extracting HGF, randomisation ' num2str(fr)])
        S.designmat = D_rand(fr).S.designmat;
        [D_rand(fr).D_fit(:).subname] = deal(D_fit.subname);
        [D_rand(fr).D_fit(:).dt] = deal(D_fit.dt);
        [out_rand(fr).T,out_rand(fr).traj,out_rand(fr).param,out_rand(fr).rt] = CORE_extract_HGF_results(D_rand(fr).D_fit,S);
    end
    
    %% stats
    [out.stats] = CORE_HGF_groupstatistics(out.T,{out.traj,S.condmean},S);
    for fr = 1:length(frand)
        disp(['calculating stats, randomisation ' num2str(fr)])
        [out_rand(fr).stats] = CORE_HGF_groupstatistics(out_rand(fr).T,{out_rand(fr).traj,S.condmean},S);
    end
    
    %% calculate randomised p values
    fields = fieldnames(out.stats.ranksum.z);
    for i = 1:length(fields)
        act_z = out.stats.ranksum.z.(fields{i});
        for rn = 1:length(out_rand)
            rand_z(rn) = [out_rand(rn).stats.ranksum.z.(fields{i})];
        end
        out.stats.ranksum.p_rand.(fields{i}) = sum(abs(rand_z) > abs(act_z)) / length(rand_z);
    end
    
    %% combined table
    if isempty(fieldnames(ranksum_all_p))
        ranksum_all_p = out.stats.ranksum.p_rand;
    else
        ranksum_all_p = [ranksum_all_p,out.stats.ranksum.p_rand];
    end
    
end
[ranksum_all_p(:).fname] = deal(fnames{:});
save(fullfile(S.path.hgf,['randomised_stats_' datestr(now,30)]),'ranksum_all_p')