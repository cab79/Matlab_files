close all
addpath('C:\Data\Matlab\export_fig')
%% figure options
fontsize = 12; % SENSITIVE TO SCREEN RESOLUTION, e.g. when using Teamviewer
save_figs=1;
savefigspath = 'C:\Users\cab79\Google Drive\1. Liverpool University\Publications\Cata anticipation 2017\Figure_images\Cue_Group_Interaction';

%% prepare SPM EEG data
S.spm_path = 'C:\Data\Catastrophising study\SPMstats\Include1\Between\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_spm';
%cluster directory name, which also specifies the constrast that will be
%plotted (i.e. the characters before the underscore)
S.clusdir='Grp_ExpAC_clusters';
%factor(s) to plot - if more than one, first factor levels will be on separate plots
%must be the same characters as used to name the factors in design_batch
S.facplot={'Exp','Grp'};
%S.facplot={'Exp'};
% clusters to plot
S.plotclus = {'c1_spm'};
S.wavetype = 'sensor'; % source or sensor?
S.wfname = 'cluster_data.mat'; %generic cluster waveform file name
S.batch = 'matlabbatch.mat'; %name of batch .mat file saved from design_batch.m and within same folder
S.subfactname = 'Subject'; %name of 'subject' factor in the SPM design
S.fact_names = {
    'Cue'
    'Group xxxxxx' % make as wide as fact names on other figure legends
    };
S.cval={ %condition labels, i.e. levels of each condition, in the same order as in the SPM design matrix. One row per factor. second column is plotting order
    {'Low','Low (prior High)','High'},[1 3]
    {'High PC','Low PC'},[1 2]
    };
S.xlimits = [-2700 -1000];% time in ms

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
St.eeglab_path = 'C:\Data\Catastrophising study\Preprocessed';
St.eventtypes = {'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'};
St.st_string='mspm12_';
St.en_string='\scond';
St.ERPsavename = 'ERP_DAT.mat';
St.overwriteEEG =0;
if ~exist(fullfile(St.eeglab_path,St.ERPsavename),'file') || St.overwriteEEG
    E=gplotprepare_eeglabdata_from_spm(St,D(1))
else
    load(fullfile(St.eeglab_path,St.ERPsavename));
end


%% set up ERP gplot
clear P
p=0;
for d = 1:length(S.plotclus);
    for sl = S.cval{1,2}
        p=p+1;
        S.selectlev = sl;
        D = gplotprepare_spmeegsensorcluster(S)
        P(p).xy = [p,1];
        P(p).x = D(d).x;
        P(p).y = D(d).y;
        P(p).cond = D(d).cond;
        P(p).condsize = []; % line/marker size
        P(p).xlinesolid = D(d).P_val;
        P(p).poly = D(d).E_val;% polygon times
        P(p).ptitle = [];%D(p).ptitle;
        P(p).fact_names = D(d).fact_names;
        P(p).colours = [0.2 0.5 1; 1 0.2 0.2]; % blue, red %CURRENTLY DOES NOT SUPPORT PLOTTING MORE THAN TWO FACTORS
        P(p).xlinedashed = [-5000 -2500 0];% vertical dashed lines to indicate events, time in ms
        P(p).timezero = -2500;% change zero latency to this time
        P(p).xaxisname = {'post-cue time (ms)'};
        P(p).yaxisname = {'amplitude (uV)'};
        P(p).plottype = 'stat_summary';
        if p~=1
            P(p).legend = 0;
        else
            P(p).legend = 1;
        end
        
        %% plot topographies
        Pt.topotype='eeglab';
        Pt.no_plot_ele=[];
        Pt.topo_subtimewin=2000;%[-2265 -2265]; % time window length. single value: plots multiple topos in consecutive windows with this length. 2 values: specifies a window. 0 = whole cluster mean.
        Pt.fontsize=fontsize;
        Pt.selectlev=S.selectlev;
        E.E_val = D(d).E_val;
        E.Fi_ind = D(d).Fi_ind;
        E.fi = D(d).fi;
        E.cond = D(d).cond;
        plot_topo(Pt,E)
        set(gcf,'color','w');
        set(gca,'fontsize', fontsize);
        set(gcf, 'Units', 'normalized', 'Position', [0.5,0.5,0.15,0.25]);
        figname = ['Topo_plot_cue_' num2str(sl)];
        if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end
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
adjustment = 1; % need to play with this to get it the same as other figures
plotheightratio = length(g)/6 / adjustment; % compared to other figures, so each subplot is the same height across figures
fig=figure;%('Position',[100 100 800 550]);
set(fig, 'Units', 'normalized', 'Position', [0,0,0.4,1*plotheightratio]);
g.set_text_options('base_size',fontsize,'legend_title_scaling',1);
g.draw();
drawnow
g=align_gplots(g,gap); % custom function
figname = 'ERP_plot';
if save_figs; export_fig(fullfile(savefigspath,[figname '.png']), '-r1200'); end

