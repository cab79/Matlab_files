close all
addpath('C:\Data\Matlab\export_fig')
%% figure options
fontsize = 12;
save_figs=1;
savefigspath = 'C:\Users\cab79\Google Drive\1. Liverpool University\Publications\Cata anticipation 2017\Figure_images\Weights';

%% Plot weight matrix scatter
clear S
S.path = 'C:\Data\Catastrophising study\SPMstats\pronto\Group Expectancy\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_GrpAvCond_gpc_ROI_perm10000';
S.fnames = {
    'weight_regress_bothExp.xlsx';
    'weight_regress_HighPCExp.xlsx';
    'weight_regress_LowPCExp.xlsx';
    };
S.ptitles = {
    'Both groups';
    'High PC group';
    'Low PC group';
    };
S.nhead = 1;
S.xvalhead = 'wimg_iv';
S.yvalhead = {'wimg_dvresid'}; % can have multiple Ys
for s = 1:length(S.fnames)
    S.fname = S.fnames{s}
    W(s) = gplotprepare_xlsxdata(S)
end
p=0;
for w = 1:length(W);
    p = p+1;
    P(p).xy = [1,p];
    P(p).x = W(w).x; % added a hack to get legend to a fixed range
    P(p).y = W(w).y;
    P(p).cond = repmat(1,length(P(p).y),1);
    P(p).condsize = repmat(1,length(P(p).y),1);
    P(p).colours = [0.4,0.4,0.4];
    P(p).ptitle = S.ptitles{p};
    P(p).xaxisname = {'Cue value weight (arbitrary units)'};
    P(p).yaxisname = {'Group weight (arbitrary units)'};
    P(p).plottype = 'geom_point_line';
    P(p).legend = 0;
end

%% draw gplot
g=gramm();
for p = 1:length(P)
    g = gplot_scatter(g,P(p));
end
fig=figure;
set(fig, 'Position', [100 100 1200 360]);
g.set_text_options('base_size',fontsize,'legend_title_scaling',1);
g.draw();
drawnow
figname = 'MVPA_weight_scatter';
if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end