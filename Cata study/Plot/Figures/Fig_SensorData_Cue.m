close all

%% prepare SPM EEG data
S.spm_path = 'C:\Data\Catastrophising study\SPMstats\Include1\Between\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_spm';
%cluster directory name, which also specifies the constrast that will be
%plotted (i.e. the characters before the underscore)
S.clusdir='ExpAC_clusters';
%factor(s) to plot - if more than one, first factor levels will be on separate plots
%must be the same characters as used to name the factors in design_batch
%P.facplot={'Exp','Grp'};
S.facplot={'Exp'};
% clusters to plot
S.plotclus = {'c1_spm','c2_spm','c3_spm','c4_spm'};
S.wavetype = 'sensor'; % source or sensor?
S.wfname = 'cluster_data.mat'; %generic cluster waveform file name
S.batch = 'matlabbatch.mat'; %name of batch .mat file saved from design_batch.m and within same folder
S.subfactname = 'Subject'; %name of 'subject' factor in the SPM design
S.fact_names = {
    'Cue'
    };
S.cval={ %condition labels, i.e. levels of each condition, in the same order as in the SPM design matrix. One row per factor. second column is plotting order
    {'Low','Low (prior High)','High'},[1 3]
    };
S.xlimits = [-2700 -1000];% time in ms
D = gplotprepare_spmeegsensorcluster(S)

%% prepare weights data for gplot
clear S
S.path = 'C:\Data\Catastrophising study\SPMstats\pronto\t-5500_1500_b-5500_-5000_m_0_1500_Grp_Exp_Subject_orig_cleaned_trialNmatch_prt_Int_gpc_ROI_noperm';
S.fname = 'weight_timeseries.xlsx';
S.nhead = 1;
S.xvalhead = 'time';
S.yvalhead = {'wimg_val','pwimg_val'};
S.condhead = {'wimg_stdev','pwimg_stdev'};
W = gplotprepare_xlsxdata(S)

%% set up ERP gplot
clear P
% Subplots 1 to 4: SPM EEG clusters 1 to 4
for p = 1:length(D);
    P(p).xy = [p,1];
    P(p).x = D(p).x;
    P(p).y = D(p).y;
    P(p).cond = D(p).cond;
    P(p).xlinesolid = D(p).P_val;
    P(p).poly = D(p).E_val;% polygon times
    P(p).ptitle = D(p).ptitle;
    P(p).fact_names = D(p).fact_names;
    P(p).colours = [0.2 0.5 1; 1 0.2 0.2]; % blue, red %CURRENTLY DOES NOT SUPPORT PLOTTING MORE THAN TWO FACTORS
    P(p).xlinedashed = [-5000 -2500 0];% vertical dashed lines to indicate events, time in ms
    P(p).timezero = -2500;% change zero latency to this time
    P(p).xaxisname = {'post-cue time (ms)'};
    P(p).yaxisname = {'amplitude (uV)'};
    P(p).plottype = 'stat_summary';
end
% Subplots 5 to 6: Weight and projection plots from MVPA
for w = 1:length(W);
    p = p+1;
    P(p).xy = [p,1];
    P(p).x = W(w).x;
    P(p).y = W(w).y;
    P(p).cond = W(w).cond;
    P(p).xlinesolid = [];
    P(p).poly = [];% polygon times
    P(p).ptitle = W(w).ptitle;
    P(p).fact_names = {'Number of SDs from mean'};
    P(p).colours = 'lch'; % blue, red %CURRENTLY DOES NOT SUPPORT PLOTTING MORE THAN TWO FACTORS
    P(p).xlinedashed = [-5000 -2500 0];% vertical dashed lines to indicate events, time in ms
    P(p).timezero = -2500;% change zero latency to this time
    P(p).xaxisname = {'post-cue time (ms)'};
    P(p).yaxisname = {'amplitude (uV)'};
    P(p).plottype = 'geom_point';
end

%% draw gplot
g=gramm();
for p = 1:length(P)
    if ~isfield(P(p),'timezero')
        P(p).timezero = [];
    end
    if ~isempty(P(p).timezero)
        P(p).xlinedashed = P(p).xlinedashed-P(p).timezero; % change zero time
        P(p).x = P(p).x-P(p).timezero;
        P(p).poly = P(p).poly-P(p).timezero;
        P(p).xlinesolid = P(p).xlinesolid-P(p).timezero;
    end
    g = gplot_timeseries(g,P(p))
end
%g.set_title('title');
g.set_text_options('base_size',4);
fig=figure;%('Position',[100 100 800 550]);
set(fig, 'Units', 'normalized', 'Position', [0,0,0.5,1]);
g.draw();

%% Prepare EEGLAB data for topos
clear S
S.eeglab_path = 'C:\Data\Catastrophising study\Preprocessed';
S.eventtypes = {'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'};
S.st_string='mspm12_';
S.en_string='\scond';
S.ERPsavename = 'ERP_DAT.mat';
if ~exist(fullfile(S.eeglab_path,S.ERPsavename),'file')
    D=gplotprepare_eeglabdata_from_spm(S,D(1))
else
    load(fullfile(S.eeglab_path,S.ERPsavename));
end

%% plot topographies
clear P
P.no_plot_ele=[];
P.topo_subtimewin=0;%[-2265 -2265]; % time window length. single value: plots multiple topos in consecutive windows with this length. 2 values: specifies a window. 0 = whole cluster mean.
plot_topo(P,D)

%% save
%psize = 18;
%g.export('file_name',[cllabel '_WF'],'export_path',fullfile(S.spm_path,S.clusdir),'file_type','svg','width',S.Nplots*psize,'height',psize,'units','centimeters');
%print(fullfile(S.spm_path,S.clusdir,[cllabel '_topo']),'-dpng');