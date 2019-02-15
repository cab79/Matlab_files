% extract, tabulate and save parameters and summary stats of trajectories

close all
clear all
dbstop if error
restoredefaultpath

% add toolbox paths
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\raincloud_plots');
addpath('C:\Data\Matlab\Violinplot-Matlab-master');
addpath('C:\Data\Matlab\cbrewer'); % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

% file info
S.path.hgf = 'C:\Data\CORE\behaviour\hgf'; 

fnames = {


% empirical priors
%'CORE_fittedparameters_percmodel12_respmodel15_fractrain1_20190106T120536_it5.mat'; %'GBM_config_BO_al4_alvar1'
%'CORE_fittedparameters_percmodel12_respmodel15_fractrain1_20190108T132736_it5.mat'; %'GBM_config_BO_al4_alvar2'
%'CORE_fittedparameters_percmodel12_respmodel15_fractrain1_20190104T200202_it7.mat'; %'GBM_config_BO_al4_alvar4'
%'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190111T203212_it2.mat'; % 'GBM_config_CORE_percmodel3_BOpriors'
%'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190111T203212_it3.mat'; % 'GBM_config_CORE_percmodel3_BOpriors'
%'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190111T203212_it4.mat'; % 'GBM_config_CORE_percmodel3_BOpriors'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it1.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it2.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it3.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it4.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it5.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it6.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it7.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it8.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it9.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it10.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it11.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it12.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it13.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it14.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it15.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it16.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it17.mat'
% 'CORE_fittedparameters_percmodel3_respmodel4_fractrain1_20190114T061324_it18.mat'

% 'CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat'
% 'D_fit_r1_it7_pm3_rm4.mat'
'D_fit_r1_it7_pm2_rm5.mat'

};

get_dt = 'CORE_fittedparameters_percmodel2_respmodel4_fractrain0_20190211T074650.mat';

ranksum_all_p = struct;
ranksum_all_z = struct;
medians_all = struct;
mva_all = struct;
for f=1:length(fnames)

    load(fullfile(S.path.hgf,'fitted',fnames{f}));
    if exist('GS','var')
        if ~isfield(D_fit,'dt')
            dtfile = load(fullfile(S.path.hgf,'fitted',get_dt));
            [D_fit(:).dt]=deal(dtfile.D_fit(:).dt);
            [D_fit(:).subname]=deal(dtfile.D_fit(:).subname);
        end
        S=GS(1);
    end
    
    % which trial-wise stats to output for traj? Averages over trials
    % without weighting by condition
    S.summary_stats = {};%{'trial_mean','trial_std','trial_absmean','trial_absstd'}; % for all traj

    % which traj to outputs to mean over conditions? First averages for
    % each condition, before creating further averages between conditions
    %S.condmean = {'PL_muhat_1','PL_dau','PL_sahat_1','PL_muhat_2','PL_sahat_2','PL_muhat_3','PL_sahat_3'};
    S.condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3','PL_psi_1','PL_psi_2','PL_psi_3','PL_sa_1','PL_sa_2','PL_sa_3','PL_sahat_1','PL_sahat_2','PL_sahat_3','PL_mu_1','PL_mu_2','PL_mu_3','PL_muhat_1','PL_muhat_2','PL_muhat_3'};
    %S.condmean = {'PL_dau','PL_da_1','PL_da_2','PL_epsi_1','PL_epsi_2','PL_epsi_3'};
    %S.condmean = {'PL_dau','PL_da_1','PL_da_2'};
    
%% for trajectory plots (means and variances on separate rows of S.condmean) 
%     S.condmean = {'PL_dau','PL_da_1','PL_da_2'; % for plotting traj: means
%             'PL_psi_1','PL_psi_2','PL_psi_3'}; var_pi = 'pi'; % precision (pi) or variance (var)?
%     S.condmean = {'PL_mu_1','PL_mu_2','PL_mu_3';
%         'PL_sa_1','PL_sa_2','PL_sa_3'}; var_pi = 'var'; % precision (pi) or variance (var)?

%%
    %S.condmean = {};

    % event numbers for each condition (1st = change, 2nd = no change)
%     S.cond.CP10_Sideaff_DC1 = [1 3]; %0.1 prob, 1 digit change, left hand
%     S.cond.CP10_Sideaff_DC3 = [2 4]; %0.1 prob, 3 digit change, left hand
%     S.cond.CP10_Sideunaff_DC1 = [5 7]; %0.1 prob, 1 digit change, right hand
%     S.cond.CP10_Sideunaff_DC3 = [6 8]; %0.1 prob, 3 digit change, right hand
%     S.cond.CP30_Sideaff_DC1 = [9 11]; %0.3 prob, 1 digit change, left hand
%     S.cond.CP30_Sideaff_DC3 = [10 12]; %0.3 prob, 3 digit change, left hand
%     S.cond.CP30_Sideunaff_DC1 = [13 15]; %0.3 prob, 1 digit change, right hand
%     S.cond.CP30_Sideunaff_DC3 = [14 16]; %0.3 prob, 3 digit change, right hand
%     S.cond.CP50_Sideaff_DC1 = [17 19]; %0.5 prob, 1 digit change, left hand
%     S.cond.CP50_Sideaff_DC3 = [18 20]; %0.5 prob, 3 digit change, left hand
%     S.cond.CP50_Sideunaff_DC1 = [21 23]; %0.5 prob, 1 digit change, right hand
%     S.cond.CP50_Sideunaff_DC3 = [22 24]; %0.5 prob, 3 digit change, right hand

%% marginal means of interest
     S.cond.condmean = [1:24]; %0.5 prob, 3 digit change, right hand
%     S.cond.odd = [1 5 9 13 17 21 2 6 10 14 18 22];  
%     S.cond.stan = [3 7 11 15 19 23 4 8 12 16 20 24];  
%      S.cond.odd_stan = {[1 5 9 13 17 21 2 6 10 14 18 22],[3 7 11 15 19 23 4 8 12 16 20 24]}; % expecting grp effects to depend on oddball vs. standard
%      S.cond.odd_CP10 = [1:2 5:6]; 
%      S.cond.odd_CP50 = [17:18 21:22]; 
%      S.cond.stan_CP10 = [3:4 7:8]; 
%      S.cond.stan_CP50 = [19:20, 23:24]; 
%      S.cond.odd_CP = {[1:2 5:6],[17:18 21:22]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_DC = {[1 5 9 13 17 21],[2 6 10 14 18 22]}; % small DC effects on odd vs. stan 
%      S.cond.stan_CP = {[3:4 7:8],[19:20, 23:24]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.stan_DC = {[3 7 11 15 19 23],[4 8 12 16 20 24]}; % small DC effects on odd vs. stan 
%      S.cond.odd_stan_CP = {[1:2 5:6],[3:4 7:8],[17:18 21:22],[19:20, 23:24]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_stan_CP10 = {[1:2 5:6],[3:4 7:8]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_stan_CP50 = {[17:18 21:22],[19:20, 23:24]}; % larger oddball effects for CP10 than CP50. Grp effects expected.
%      S.cond.odd_stan_DC = {[1 5 9 13 17 21],[3 7 11 15 19 23],[2 6 10 14 18 22],[4 8 12 16 20 24]}; % small DC effects on odd vs. stan 

%% marginal means for violin plots
%     S.cond.CP10_Odd_DC1 = [1 5];
%     S.cond.CP10_Odd_DC3 = [2 6];
%     S.cond.CP30_Odd_DC1 = [9 13];
%     S.cond.CP30_Odd_DC3 = [10 14];
%     S.cond.CP50_Odd_DC1 = [17 21];
%     S.cond.CP50_Odd_DC3 = [18 22];
%     S.colour_code = [1 2 1 2 1 2];
%     S.x_pos = [1.2 1.8 3.2 3.8 5.2 5.8]; % x position of each violin
%     S.abs = 0;

%% each cell of design: for outputting to excel for SPSS, or for plotting
%       S.cond.CP10_SideA_Odd_DC1 = 1;
%       S.cond.CP10_SideA_Odd_DC3 = 2;
%       S.cond.CP10_SideA_Stan_DC1 = 3;
%       S.cond.CP10_SideA_Stan_DC3 = 4;
%       S.cond.CP10_SideU_Odd_DC1 = 5;
%       S.cond.CP10_SideU_Odd_DC3 = 6;
%       S.cond.CP10_SideU_Stan_DC1 = 7;
%       S.cond.CP10_SideU_Stan_DC3 = 8;
%       S.cond.CP30_SideA_Odd_DC1 = 9;
%       S.cond.CP30_SideA_Odd_DC3 = 10;
%       S.cond.CP30_SideA_Stan_DC1 = 11;
%       S.cond.CP30_SideA_Stan_DC3 = 12;
%       S.cond.CP30_SideU_Odd_DC1 = 13;
%       S.cond.CP30_SideU_Odd_DC3 = 14;
%       S.cond.CP30_SideU_Stan_DC1 = 15;
%       S.cond.CP30_SideU_Stan_DC3 = 16;
%       S.cond.CP50_SideA_Odd_DC1 = 17;
%       S.cond.CP50_SideA_Odd_DC3 = 18;
%       S.cond.CP50_SideA_Stan_DC1 = 19;
%       S.cond.CP50_SideA_Stan_DC3 = 20;
%       S.cond.CP50_SideU_Odd_DC1 = 21;
%       S.cond.CP50_SideU_Odd_DC3 = 22;
%       S.cond.CP50_SideU_Stan_DC1 = 23;
%       S.cond.CP50_SideU_Stan_DC3 = 24;
%     S.colour_code = repmat([1 2],1,12);
%     S.x_pos = [1:8 13:20 25:32]; % x position of each violin
%     S.abs = 0;

%%     
%     S.cond.condmean_odd_stan = {[1 2 5 6 9 10 13 14 17 18 21 22],[3 4 7 8 11 12 15 16 19 20 23 24]}; %0.5 prob, 3 digit change, right hand
%     S.cond.DC1_odd_stan = {[1 5 9 13 17 21],[3 7 11 15 19 23]}; %0.5 prob, 3 digit change, right hand
%     S.cond.DC3_odd_stan = {[2 6 10 14 18 22],[4 8 12 16 20 24]}; %0.5 prob, 3 digit change, right hand
%     S.cond.CP10_odd_stan = {[1 2 5 6],[3 4 7 8]}; %0.5 prob, 3 digit change, right hand
%     S.cond.CP30_odd_stan = {[9 10 13 14],[11 12 15 16]}; %0.5 prob, 3 digit change, right hand
%     S.cond.CP50_odd_stan = {[17 18 21 22],[19 20 23 24]}; %0.5 prob, 3 digit change, right hand
%     S.cond.sideaff_odd_stan = {[1 2 9 10 17 18],[3 4 11 12 19 20]}; %0.5 prob, 3 digit change, right hand
%     S.cond.sideunaff_odd_stan = {[5 6 13 14 21 22],[7 8 15 16 23 24]}; %0.5 prob, 3 digit change, right hand
    
    % 
    % S.cond.CP10_Sideaff_DC1_ch = [1]; %0.1 prob, 1 digit change, left hand
    % S.cond.CP10_Sideaff_DC3_ch = [2]; %0.1 prob, 3 digit change, left hand
    % S.cond.CP10_Sideunaff_DC1_ch = [5]; %0.1 prob, 1 digit change, right hand
    % S.cond.CP10_Sideunaff_DC3_ch = [6]; %0.1 prob, 3 digit change, right hand
    % S.cond.CP30_Sideaff_DC1_ch = [9]; %0.3 prob, 1 digit change, left hand
    % S.cond.CP30_Sideaff_DC3_ch = [10]; %0.3 prob, 3 digit change, left hand
    % S.cond.CP30_Sideunaff_DC1_ch = [13]; %0.3 prob, 1 digit change, right hand
    % S.cond.CP30_Sideunaff_DC3_ch = [14]; %0.3 prob, 3 digit change, right hand
    % S.cond.CP50_Sideaff_DC1_ch = [17]; %0.5 prob, 1 digit change, left hand
    % S.cond.CP50_Sideaff_DC3_ch = [18]; %0.5 prob, 3 digit change, left hand
    % S.cond.CP50_Sideunaff_DC1_ch = [21]; %0.5 prob, 1 digit change, right hand
    % S.cond.CP50_Sideunaff_DC3_ch = [22]; %0.5 prob, 3 digit change, right hand
    % 
    % S.cond.CP10_Sideaff_DC1_noch = [3]; %0.1 prob, 1 digit change, left hand
    % S.cond.CP10_Sideaff_DC3_noch = [4]; %0.1 prob, 3 digit change, left hand
    % S.cond.CP10_Sideunaff_DC1_noch = [7]; %0.1 prob, 1 digit change, right hand
    % S.cond.CP10_Sideunaff_DC3_noch = [8]; %0.1 prob, 3 digit change, right hand
    % S.cond.CP30_Sideaff_DC1_noch = [11]; %0.3 prob, 1 digit change, left hand
    % S.cond.CP30_Sideaff_DC3_noch = [12]; %0.3 prob, 3 digit change, left hand
    % S.cond.CP30_Sideunaff_DC1_noch = [15]; %0.3 prob, 1 digit change, right hand
    % S.cond.CP30_Sideunaff_DC3_noch = [16]; %0.3 prob, 3 digit change, right hand
    % S.cond.CP50_Sideaff_DC1_noch = [19]; %0.5 prob, 1 digit change, left hand
    % S.cond.CP50_Sideaff_DC3_noch = [20]; %0.5 prob, 3 digit change, left hand
    % S.cond.CP50_Sideunaff_DC1_noch = [23]; %0.5 prob, 1 digit change, right hand
    % S.cond.CP50_Sideunaff_DC3_noch = [24]; %0.5 prob, 3 digit change, right hand


    %% run stats
    if ~exist('D_fit','var') && exist('D_sim','var')
        D_fit=D_sim;
        D_fit.HGF.fit = D_fit.HGF.sim;
    end
    [out.T,out.traj,out.param,out.rt] = CORE_extract_HGF_results(D_fit,S);
    
    %% save table
    if 0
        sname = strrep(fullfile(S.path.hgf,'fitted',fnames{f}),'.mat','_table.xlsx');
        writetable(out.T,sname)
    end
    
    if length(D_fit)>1
        S.lda = {'traj_conds'};%{'para','traj_trials','traj_conds'};
        S.nperm=0;
        %[out.stats] = CORE_HGF_groupstatistics(out.T,{out.traj},S);
        [out.stats] = CORE_HGF_groupstatistics(out.T,{out.traj,S.condmean},S);
        %[out.stats] = CORE_HGF_groupstatistics(out.T,{},S);
    end
    
    %% combined table
    if isempty(fieldnames(ranksum_all_p))
        ranksum_all_p = out.stats.ranksum.p;
    else
        ranksum_all_p = [ranksum_all_p,out.stats.ranksum.p];
    end
    
    if isempty(fieldnames(ranksum_all_z))
        ranksum_all_z = out.stats.ranksum.z;
    else
        ranksum_all_z = [ranksum_all_z,out.stats.ranksum.z];
    end
    
    %% combined table
    try
    if isempty(fieldnames(mva_all))
        mva_all = out.stats.mvc.traj.error;
    else
        mva_all = [mva_all,out.stats.mvc.traj.error];
    end
    end
    
    %% group medians
    [grps,~,grpind] = unique(out.T.groups);
    hdr = out.T.Properties.VariableNames;
    out.median = table;
    for i = 1:size(out.T,2)
        if isnumeric(out.T.(hdr{i})(1))
            for g = 1:length(grps)
                out.median.(hdr{i})(g) = median(out.T.(hdr{i})(grpind==g));
            end
        end
    end
    %% combined table
    if isempty(fieldnames(medians_all))
        medians_all = table2struct(out.median);
    else
        medians_all = [medians_all,table2struct(out.median)];
    end

    %% save
    if 0
        save(fullfile(S.path.hgf,strrep(fnames{f},'fittedparameters','groupstatistics')), 'out');
        sname=strrep(strrep(fname,'fittedparameters','statstable'),'.mat','.xlsx');
        writetable(out.T,fullfile(S.path.hgf,sname));
    end

    %% plot group differences for selected variables
    if 1
        close all
        %vs={'PL_dau_odd','PL_dau_stan','PL_dau_odd_stan','PL_dau_odd_CP','PL_dau_stan_CP','PL_dau_odd_stan_CP','PL_epsi_1_odd','PL_epsi_2_odd','PL_epsi_3_odd','PL_epsi_1_stan','PL_epsi_2_stan','PL_epsi_3_stan','PL_psi_1_condmean','PL_psi_2_condmean','PL_psi_3_condmean','PL_psi_1_odd','PL_psi_2_odd','PL_psi_3_odd','PL_psi_1_stan','PL_psi_2_stan','PL_psi_3_stan','PL_epsi_1_odd_CP','PL_epsi_1_stan_CP','PL_epsi_1_odd_stan_CP','PL_dau_odd_DC','PL_dau_stan_DC','PL_dau_odd_stan_DC','PL_da_1_odd','PL_da_1_odd_CP','PL_da_1_odd_DC','PL_epsi_2_odd_CP','PL_epsi_2_stan_CP','PL_epsi_2_odd_stan_CP','PL_da_2_odd','PL_da_2_odd_CP','PL_da_2_odd_DC','PL_epsi_3_odd_CP','PL_epsi_3_stan_CP','PL_epsi_3_odd_stan_CP'};
        %vs={'like_al0_1','like_al0_2','like_al0_3','like_al0_4','PL_om_2','PL_om_3','PL_dau_mean','PL_da_1_mean','PL_da_2_mean','PL_epsi_1_mean','PL_epsi_2_mean','PL_epsi_3_mean','PL_psi_1_mean','PL_psi_2_mean','PL_psi_3_mean','PL_dau_odd','PL_dau_stan','PL_dau_odd_stan','PL_dau_odd_CP','PL_dau_stan_CP','PL_dau_odd_stan_CP','PL_dau_odd_DC','PL_dau_stan_DC','PL_dau_odd_stan_DC','PL_da_1_odd','PL_da_1_odd_CP','PL_da_1_odd_DC','PL_da_2_odd','PL_da_2_odd_CP','PL_da_2_odd_DC'};
        %vs={'like_al0','PL_om_2','PL_om_3','PL_sa_2_mean','PL_ud_1_mean','PL_ud_2_mean','PL_dau_mean','PL_epsi_2_mean','PL_epsi_3_mean'};
        vs={'like_al0_1','like_al0_2','like_al0_3','PL_om_2','PL_om_3','PL_dau_condmean','PL_da_1_condmean','PL_da_2_condmean','PL_epsi_1_condmean','PL_epsi_2_condmean','PL_epsi_3_condmean','PL_mu_1_condmean','PL_mu_2_condmean','PL_mu_3_condmean','PL_sa_1_condmean','PL_sa_2_condmean','PL_sa_3_condmean'};
        titles = {'alpha','alphaSide','alphaDC','omega','theta','dau','da 1','da 2','epsi 1','epsi 2','epsi 3','mu 1','mu 2','mu 3','sa 1','sa 2','sa 3'};
        sp1 = 6; % rows
        sp2 = 3; % columns
        vn=0;figure
        for v=1:length(vs)
            vn=vn+1;
            subplot(sp1,sp2,vn)
            for g = 1:length(grps)
                % average over cell elements of vs
                if iscell(vs{v})
                    for vii = 1:length(vs{v})
                        X{g}(:,vii)=out.T.(vs{v}{vii})(grpind==g);
                    end
                    X{g} = mean(X{g},2);
                else
                    X{g}=out.T.(vs{v})(grpind==g);
                end
            end
            [cb] = cbrewer('qual', 'Set1', g, 'pchip');
            h = n_rainclouds(X,cb);
            if exist('titles','var')
                title(titles{v})
            else
                title([vs{v} ', p = ' num2str(out.stats.ranksum.p.(vs{v}))])
            end
            set(gca,'YTick',[],'fontsize',10,'ycolor',[1 1 1])
            box off
        end
        legend(grps)
    end
    
    %% HGF trajectories
    if 0
        close all
        subplot_on = 1;
        plot_variance_type = 'mean_of_traj'; %{'mean_of_traj','var_over_subs'};
        
        % group indices
        for g = 1:length(grps)
            % Set up display
            scrsz = get(0,'screenSize');
            outerpos = [0.2*scrsz(3),0.2*scrsz(4),0.8*scrsz(3),0.8*scrsz(4)];
            figure(...
                'OuterPosition', outerpos,...
                'Name', ['HGF trajectories, ' grps{g}]);

            vs=S.condmean(1,:);
            [cb] = cbrewer('qual', 'Set1', length(vs), 'pchip');

            leg={''};
            for v=1:length(vs)

                if subplot_on
                    ax(g,v)=subplot(length(vs),1,length(vs)-v+1);
                else
                    hold on
                end

                % get data
                dat = cat(2,out.traj(:).(vs{v}));
                mn = mean(dat(:,grpind==g),2)';
                ns = length(grpind==g);
                ts=1:length(mn);

                if strcmp(plot_variance_type,'var_over_subs')
                    % get spread
                    sd = std(dat(:,grpind==g),[],2)';
                    sem = sd/sqrt(ns);              % Standard Error
                    tscore = tinv(0.05,ns-1);      % T-Score
                    CI = tscore*sem;                % Confidence Intervals
                    % choose spread type
                    spr = sd/2;
                elseif strcmp(plot_variance_type,'mean_of_traj')
                    vs_var=S.condmean{2,v};
                    % get data
                    dat_var = cat(2,out.traj(:).(vs_var));
                    dat_var = dat_var(:,grpind==g);
                    if strcmp(var_pi,'pi')
                        dat_var = 1./dat_var;
                    end
                    spr = mean(dat_var,2)'/2;
                end

                upper = mn+spr;
                lower = mn-spr;

                % plot spread
                fill([ts, fliplr(ts)], [(upper), fliplr((lower))], ...
                cb(v,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.15);
                hold on

                % plot mean
                ln(v) = plot(ts,mn,'Color',cb(v,:));
                axval = axis;
                lim(g,v,:)=[axval(1:2) min(lower) max(upper)];
                if subplot_on
                    ylabel(vs{v})
                    if v==1
                        xlabel('Trial number')
                    else
                        set(gca, 'XTickLabel', [])
                    end
                end
                hold off
                set(gca,'FontSize',18)
                %tightfig
            end
            % plot inputs/design
            %plot(ts, D_fit(1).HGF.u(:,1), '.', 'Color', [0 0.6 0]); % inputs
            if ~subplot_on
                legend(ln,vs)
                ylabel('Prediction error')
                xlabel('Trial number')
            end
        end
        % equate axes between figures
        maxlim(:,[2,4]) = squeeze(max(lim(:,:,[2 4]),[],1));
        maxlim(:,[1,3]) = squeeze(min(lim(:,:,[1 3]),[],1));
        for g = 1:length(grps)
            for v=1:length(vs)
                axis(ax(g,v),maxlim(v,:));
            end
        end
    end
    


    %% condition means of HGF trajectories
    if 0
        vs=S.condmean;
        ucol = unique(S.colour_code);
        [cb] = cbrewer('qual', 'Set1', length(ucol), 'pchip');
        for v=1:length(vs)

            % get column indices
            colnames = out.T.Properties.VariableNames;
            vss{v} = colnames(~cellfun(@isempty,strfind(colnames,vs{v})));
            
            % get data
            dat = out.T(:,vss{v});
            xlab = dat.Properties.VariableNames;
            if S.abs
                plotdat = abs(dat{:,:});
            else
                plotdat = dat{:,:};
            end

            % plot
            figure('Name','All groups')
            %bar(xlab,dat{1,:})
            violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[0 0 0])
            title(vs{v})
            
            % separate plots by group indices
            for g = 1:length(grps)

                % get data
                dat = out.T(grpind==g,vss{v});
                xlab = dat.Properties.VariableNames;
                if S.abs
                    plotdat = abs(dat{:,:});
                else
                    plotdat = dat{:,:};
                end

                % plot
                figure('Name',grps{g})
                %bar(xlab,dat{1,:})
                violinplot(plotdat,xlab,S.x_pos,'ViolinColor',cb(S.colour_code,:),'ViolinAlpha',0.6,'EdgeColor',[0.8 0.8 0.8],'MedianColor',[0 0 0])
                title(vs{v})
            end
        end
    end
end
try [mva_all(:).fname] = deal(fnames{:});end
[ranksum_all_p(:).fname] = deal(fnames{:});

% flag significance
fd=fieldnames(ranksum_all_p);
ii = length(ranksum_all_p)+1;
for fn = 1:length(fd)
    ranksum_all_p(ii).(fd{fn}) = sum([ranksum_all_p(1:ii-1).(fd{fn})]<0.05)>length(ranksum_all_p)/10;
end
