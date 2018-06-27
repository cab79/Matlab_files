clear all
close all
addpath('C:\Data\Matlab\export_fig')
%% figure options
fontsize = 12; % SENSITIVE TO SCREEN RESOLUTION, e.g. when using Teamviewer
save_figs=1;
%savefigspath = 'C:\Data\CORE\EEG\ana\spm\SPMstats\t-200_899_b-200_0_m_0_600_Grp_Odd_DC_Subject_2_merged_cleaned_spm\Grp_clusters';

%% prepare SPM EEG data
S.spm_path = 'C:\Data\CORE\EEG\ana\spm\SPMstats\t-200_899_b-200_0_m_0_600_Grp_Odd_DC_Subject_2_merged_cleaned_spm';
%cluster directory name, which also specifies the constrast that will be
%plotted (i.e. the characters before the underscore)
S.clusdir='Odd_clusters';
%factor(s) to plot - if more than one, first factor levels will be on separate plots
%must be the same characters as used to name the factors in design_batch
S.facplot={'Odd'};
%S.facplot={'Exp'};
% clusters to plot
S.plotclus = {'c1_spm','c2_spm','c3_spm'};
S.plotclus_sep = 1; % separate plots for each cluster
S.wavetype = 'sensor'; % source or sensor?
S.wfname = 'cluster_data.mat'; %generic cluster waveform file name
S.batch = 'matlabbatch.mat'; %name of batch .mat file saved from design_batch.m and within same folder
S.subfactname = 'Subject'; %name of 'subject' factor in the SPM design
S.fact_names = {
    %'Change Probability';
    %'Group';
    'Oddball effect';
    %'Digit Change';
    %'Side';
    };
S.cval={ %condition labels, i.e. levels of each condition, in the same order as in the SPM design matrix. One row per factor. second column is plotting order
    %{'10%','30%','50%'};
    %{'CRPS','HC'},[1 2];
    {'Oddball','Standard'},[1 2];
    %{'DC1','DC3'};
    %{'Affected','Unaffected'}
    };
S.xlimits = [-200 800];% time in ms
D = gplotprepare_spmeegsensorcluster(S)

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
St.st_string='mspm12_';
St.en_string='\scond';
St.ERPsavename = 'ERP_DAT.mat';
St.plot_diff_wave = 1;
St.use_flipped=1;
St.flipcond = [5:8];
St.overwriteEEG =0;
if ~exist(fullfile(St.eeglab_path,St.ERPsavename),'file') || St.overwriteEEG
    E=gplotprepare_eeglabdata_from_spm(St,D(1))
else
    load(fullfile(St.eeglab_path,St.ERPsavename));
    for fn = fieldnames(D)'
       E.(fn{1}) = D.(fn{1});
    end
end


%% set up ERP gplot
clear P
% Subplots 1 to 4: SPM EEG clusters 1 to 4
for p = 1:length(D);
    P(p).xy = [p,1];
    P(p).x = D(p).x;
    P(p).y = D(p).y;
    P(p).cond = D(p).cond;
    P(p).condsize = []; % line/marker size
    P(p).xlinesolid = D(p).P_val;
    P(p).poly = D(p).E_val;% polygon times
    P(p).ptitle = [];%D(p).ptitle;
    P(p).fact_names = D(p).fact_names;
    P(p).colours = [0.2 0.5 1; 1 0.2 0.2]; % blue, red %CURRENTLY DOES NOT SUPPORT PLOTTING MORE THAN TWO FACTORS
    P(p).xlinedashed = [0];% vertical dashed lines to indicate events, time in ms
    P(p).timezero = 0;% change zero latency to this time
    P(p).xaxisname = {'post-cue time (ms)'};
    P(p).yaxisname = {'amplitude (uV)'};
    P(p).plottype = 'stat_summary';
    if p~=1
        P(p).legend = 0;
    end
end
%% draw gplot
g=gramm();
for p = 1:length(P)
    if ~isfield(P(p),'timezero'); 
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

%g.set_title('title');
gap = 1.2; % 1.3 = 30% gap between plots
adjustment = 1.05; % need to play with this to get it the same as other figures
plotheightratio = length(g)/6 / adjustment; % compared to other figures, so each subplot is the same height across figures
fig=figure;%('Position',[100 100 800 550]);
set(fig, 'Units', 'normalized', 'Position', [0,0,0.4,1*plotheightratio]);
g.set_text_options('base_size',fontsize,'legend_title_scaling',1);
g.draw();
drawnow
g=align_gplots(g,gap); % custom function
figname = 'EEG_plot';
if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end

%% plot topographies
clear P
P.topotype='eeglab';
P.no_plot_ele=[];
P.topo_subtimewin=2000;%[-2265 -2265]; % time window length. single value: plots multiple topos in consecutive windows with this length. 2 values: specifies a window. 0 = whole cluster mean.
P.fontsize=fontsize;
P.cval = 0;
for d = 1:length(D)
    E.E_val = D(d).E_val
    plot_topo(P,E)
    set(gcf,'color','w');
    set(gca,'fontsize', fontsize);
    set(gcf, 'Units', 'normalized', 'Position', [0.5,0.5,0.15,0.25]);
    figname = [S.plotclus{d} '_topo_plot_' num2str(E.E_val(1)) '-' num2str(E.E_val(end))];
    if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end
end


