% Searchlight analysis
% RECOMMEND TO USE WITH NO MORE THAN 2-FOLD CROSS-VALIDATION FOR SPEED
clear all
Din = 'C:\Data\Catastrophising study\SPMstats\pronto\t-3000_-2_b-3000_-2500_m_-2500_-1000_Exp_orig_cleaned_SPNall_prt_Exp_gpc_ROI_lobo';
Din_sub = 'sub*';
%opt.searchtype = 'spacetime'; opt.R = 5;
opt.searchtype = 'time'; opt.R = [inf inf 100];
opt.parallel = 1;
opt.i_model = 1; % Model index to use
opt.loadF = 1; % load all features
opt.savImg = 0;
opt.plot = 0;
opt.permStat = 0; % permutations

%% RUN
subdir = dir(fullfile(Din,Din_sub));
dirFlags = [subdir.isdir];
subdir = subdir(dirFlags);
for s = 1:length(subdir)
    Sub_in = fullfile(Din,subdir(s).name,'PRT.mat');
    [SLres{s},Pout{s},XYZ] = crc_parSL(Sub_in,opt);
end

%% Group info
load(fullfile(Din,'sub_info.mat'))
grp=[];
for g = 1:length(SubInd)
    grp = [grp g*ones(1,length(SubInd{g}))];
end

%% PLOT
figure
hold on
cols={'r','b'};
for s = 1:length(SLres)
    acc(:,s) = [SLres{s}(:).acc];
    col=cols{grp(s)};
    plot(XYZ(3,:),acc(1:size(XYZ,2),s),col); 
    acc_lb(:,s) = [SLres{s}(:).acc_lb];
    acc_ub(:,s) = [SLres{s}(:).acc_ub];
    %plot(XYZ(3,:),acc_lb(:,s),'--')
    %plot(XYZ(3,:),acc_ub(:,s),'--')
end
hold off
for g = 1:length(SubInd)
    acc_av(:,g) = mean(acc(:,grp==g),2);
    acc_lb_av(:,g) = mean(acc_lb(:,grp==g),2);
    acc_ub_av(:,g) = mean(acc_ub(:,grp==g),2);
end
figure
hold on
plot(XYZ(3,:),acc_av(1:size(XYZ,2),:)); 
plot(XYZ(3,:),acc_lb_av(1:size(XYZ,2),:),'--')
plot(XYZ(3,:),acc_ub_av(1:size(XYZ,2),:),'--')
hold off


%% STATS
addpath(genpath('C:\Data\Matlab\TFCE'));
% independent sample ttest, two-sided
S.analysis = 'independent';
S.tails = 2;
S.g1 = find(grp==1);
S.g2 = find(grp==2);
S.imgs(:,1,1,:) = acc(:,S.g1);
S.imgs2(:,1,1,:) = acc(:,S.g2);
S.covariate = [];
S.nuisance = [];
S.H = [];
S.E = [];
S.C = [];
S.dh = [];
S.parworkers = [];
S.nperm = 5000;
[S.pcorr_pos,S.pcorr_neg] = matlab_tfce(S.analysis,S.tails,S.imgs,S.imgs2,S.covariate,S.nperm,S.H,S.E,S.C,S.dh,S.parworkers,S.nuisance);
%fp = sum(pcorr_pos(:)<.05)+sum(pcorr_neg(:)<.05)
save(fullfile(Din,'TFCE_stats.mat'),'S');
find(S.pcorr_pos<0.05)
find(S.pcorr_neg<0.05)
