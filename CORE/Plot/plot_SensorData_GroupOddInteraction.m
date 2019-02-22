clear all
close all
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
addpath('C:\Data\Matlab\export_fig')
addpath('C:\Data\CORE\eeg')
addpath(genpath('C:\Data\Matlab\gramm-master'))
%% figure options
fontsize = 12; % SENSITIVE TO SCREEN RESOLUTION, e.g. when using Teamviewer
save_figs=1;
%savefigspath = 'C:\Data\CORE\EEG\ana\spm\SPMstats\t-200_899_b-200_0_m_0_600_Grp_Odd_DC_Subject_2_merged_cleaned_spm\Grp_Odd_clusters';

%% prepare SPM EEG data
%S.spm_path = 'C:\Data\CORE\EEG\ana\spm\SPMstats\t-200_899_b-200_0_m_0_600_Grp_Odd_DC_Subject_2_merged_cleaned_spm';
S.spm_path = 'C:\Data\CORE\eeg\ana\spm\SPMstats\sensor\t-200_899_b-200_0_m_0_800_Side_Grp_Odd_Subject_Age_2_merged_cleaned_spm';
%cluster directory name, which also specifies the constrast that will be
%plotted (i.e. the characters before the underscore)
S.clusdir='Grp_Odd_clusters';
%factor(s) to plot - if more than one, first factor levels will be on separate plots
%must be the same characters as used to name the factors in design_batch
S.facplot={'Grp','Odd'};
%S.facplot={'Exp'};
% clusters to plot
S.plotclus = {'c4_spm','c3_spm','c1_spm','c5_spm','c2_spm'};
S.plotclus_sep = 1; % separate plots for each cluster
S.subplot_max = 4;
S.wavetype = 'sensor'; % source or sensor?
S.wfname = 'cluster_data.mat'; %generic cluster waveform file name
S.batch = 'matlabbatch.mat'; %name of batch .mat file saved from design_batch.m and within same folder
S.subfactname = 'Subject'; %name of 'subject' factor in the SPM design
S.fact_names = {
    %'Change Probability';
    'Group';
    'Oddball effect';
    %'Digit Change';
    %'Side';
    };
S.cval={ %condition labels, i.e. levels of each condition, in the same order as in the SPM design matrix. One row per factor. second column is plotting order
    %{'10%','30%','50%'};
    {'CRPS','HC'},[1 2];
    {'Oddball','Standard'},[1 2];
    %{'DC1','DC3'};
    %{'Affected','Unaffected'}
    };
S.xlimits = [-200 800];% time in ms
D = gplotprepare_spmeegsensorcluster(S);

%% prepare weights data for gplot
%Sw.path = 'C:\Data\Catastrophising study\SPMstats\pronto\Group\Main effect\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_prt_GrpAvCond_gpc_ROI_noperm';
%Sw.fname = 'weight_timeseries.xlsx';
%Sw.nhead = 1;
%Sw.xvalhead = 'time';
%Sw.yvalhead = {'wimg_val','pwimg_val'};
%Sw.condhead = {'wimg_stdev','pwimg_stdev'};
%Sw.colormaprange = [-3:3]; % should include all cond values 
%W = gplotprepare_xlsxdata(Sw)


%% Prepare EEGLAB data for topos
savefigspath = fullfile(S.spm_path,S.clusdir);
St.eeglab_path = 'C:\Data\CORE\EEG\ana\prep\cleaned\part2';
St.eventtypes = [1:8];%{'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'};
St.mark = {
    [1 9 17]; % left, Odd, DC1
    [2 10 18]; % left, Odd, DC3
    [3 11 19]; % left, stan, DC1
    [4 12 20]; % left, stan, DC3
    [5 13 21]; % right, Odd
    [6 14 22]; % right, Odd
    [7 15 23]; % right, stan
    [8 16 24]; % right, Stan
};
St.st_string='mspm12_flip_CPavg_';
St.en_string='\scond';
St.ERPsavename = 'ERP_DAT_CPavg.mat';
St.plot_diff_wave = 1;
St.use_flipped=1;
St.flipcond = [5:8 13:16 21:24]; % EEGLAB markers
St.overwriteEEG =0;
if ~exist(fullfile(St.eeglab_path,St.ERPsavename),'file') || St.overwriteEEG
    E=gplotprepare_eeglabdata_from_spm(St,D(1));
else
    load(fullfile(St.eeglab_path,St.ERPsavename));
    for fn = fieldnames(D)'
       E.(fn{1}) = D.(fn{1});
    end
end



%% set up ERP gplot
clear P
p=0;
for d = 1:length(S.plotclus)
    for sl = S.cval{1,2} % for each group
        p=p+1;
        S.selectlev = sl;
        D = gplotprepare_spmeegsensorcluster(S);
        P(p).xy = [p,1];
        P(p).x = D(d).x;
        P(p).y = D(d).y;
        P(p).cond = D(d).cond;
        P(p).condsize = []; % line/marker size
        P(p).xlinesolid = D(d).P_val;
        P(p).poly = D(d).E_val;% polygon times
        P(p).ptitle = [];%D(p).ptitle;
        P(p).fact_names = D(d).fact_names;
        %P(p).colours = [0.2 0.5 1; 1 0.2 0.2]; % blue, red 
        P(p).colours = [0.75 0 0; 0.25 0 0]; % purple, orange, green 
        P(p).color_order = -1; % order of plotting: -1 or 1
        P(p).xlinedashed = [0];% vertical dashed lines to indicate events, time in ms
        P(p).timezero = 0;% change zero latency to this time
        P(p).xaxisname = {'post-cue time (ms)'};
        P(p).yaxisname = {'amplitude (uV)'};
        P(p).plottype = 'stat_summary';
        %if p~=1
        %    P(p).legend = 0;
        %else
            P(p).legend = 1;
        %end
        
        %% plot topographies
        Pt.topotype='eeglab';
        Pt.no_plot_ele=[];
        Pt.topo_subtimewin=2000;%[-2265 -2265]; % time window length. single value: plots multiple topos in consecutive windows with this length. 2 values: specifies a window. 0 = whole cluster mean.
        Pt.fontsize=fontsize;
        Pt.selectlev=S.selectlev;
        Pt.cval = S.cval(2,:); % conditions
        Pt.sub_order = -1;
        E.E_val = D(d).E_val;
        E.Fi_ind = D(d).Fi_ind;
        E.fi = D(d).fi;
        E.cond = D(d).cond;
        plot_topo(Pt,E)
        set(gcf,'color','w');
        set(gca,'fontsize', fontsize);
        set(gcf, 'Units', 'normalized', 'Position', [0.5,0.5,0.15,0.25]);
        figname = [S.plotclus{d} '_topo_plot_oddcond_' num2str(sl)];
        if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end
    end
    
end

%% draw gplot
g=gramm();
for p = 1:length(P)
    if ~isfield(P(p),'timezero')
        P(p).timezero = [];
    end
    if ~isempty(P(p).timezero)
        P(p).xlinedashed = P(p).xlinedashed-P(p).timezero; % change zero time
        P(p).xlinesolid = P(p).xlinesolid-P(p).timezero;
        P(p).x = P(p).x-P(p).timezero;
        P(p).poly = P(p).poly-P(p).timezero;
    end
    g = gplot_timeseries(g,P(p));
end
n_fig = ceil(length(g)/S.subplot_max);
G=g;
for nf = 1:n_fig
    ind = (nf-1)*S.subplot_max + 1 : min((nf-1)*S.subplot_max + S.subplot_max, length(G));
    g = G(ind);
    %g.set_title('title');
    gap = 1.2; % 1.3 = 30% gap between plots
    adjustment = 1; % need to play with this to get it the same as other figures
    plotheightratio = length(g)/S.subplot_max / adjustment; % compared to other figures, so each subplot is the same height across figures
    fig=figure;%('Position',[100 100 800 550]);
    set(fig, 'Units', 'normalized', 'Position', [0,0,0.4,1*plotheightratio]);
    g.set_text_options('base_size',fontsize,'legend_title_scaling',1);
    g.draw();
    drawnow
    g=align_gplots(g,gap); % custom function
    %figname = [S.plotclus{p} '_ERP_plot'];
    figname = 'ERP_plot';
    if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end
end



