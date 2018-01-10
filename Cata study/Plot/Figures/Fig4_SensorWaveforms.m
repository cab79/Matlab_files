%% plot ERP waveforms from clusters using the amazing Gramm toolbox!

clear all
close all
%% generic directories for all analyses for this study
%-------------------------------------------------------------
% directory in which SPM is saved
%P.spm_path = 'C:\Data\Catastrophising study\SPMstats\Source\1_grp\_Time_Int_Exp_Subject_spm_t416_478';
%P.wavetype = 'source'; % source or sensor?
P.spm_path = 'C:\Data\Catastrophising study\SPMstats\Include1\Between\t-3000_-2_b-3000_-2500_m_-2500_0_Grp_Exp_Subject_orig_cleaned_SPNall_spm';
P.wavetype = 'sensor'; % source or sensor?

%% specific directory and file information for this analysis
%-------------------------------------------------------------
%generic cluster waveform file name
P.wfname = 'cluster_data.mat';
%name of batch .mat file saved from design_batch.m and within same folder
%as SPM.mat
P.batch = 'matlabbatch.mat';
%name of 'subject' factor in the SPM design
P.subfactname = 'Subject';

%% plot options
%-------------------------------------------------------------
%CURRENTLY DOES NOT SUPPORT PLOTTING MORE THAN TWO FACTORS
P.colours = [0.2 0.5 1; 1 0.2 0.2]; % blue, red
% axis names - can differ according the experiment
P.xaxisname = {'peri-stimulus time (ms)'};
P.yaxisname = {'amplitude (arbitrary units)'};
% vertical dashed lines to indicate events, time in ms
P.xlimits = [];
% vertical dashed lines to indicate events, time in ms
P.xlines = [-5000 -2500 0];
% polygon times (or leave blank to extract from clustertable)
P.poly = [];

% plot topographies
P.plot_topo=1;
P.eventtypes = {'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'};
P.eeglab_path = 'C:\Data\Catastrophising study\Preprocessed';
P.no_plot_ele=[];
P.topo_subtimewin=[-2265 -2265]; % time window length. single value: plots multiple topos in consecutive windows with this length. 2 values: specifies a window. 0 = whole cluster.

% save optons
P.save_waveforms=0;
P.save_topo=1;

%% figures to plot
P.plotnum = 3;

switch P.plotnum
    case 1
        %cluster directory name, which also specifies the constrast that will be
        %plotted (i.e. the characters before the underscore)
        P.clusdir='ExpAB_clusters';
        %factor(s) to plot - if more than one, first factor levels will be on separate plots
        %must be the same characters as used to name the factors in design_batch
        %P.facplot={'Int','Exp'};
        %P.facplot={'Exp','Grp'};
        %P.facplot={'Att','Exp'};
        P.facplot={'Exp'};
        %Full factor names to use for labelling the plots
        P.fact_names = {
            %'Attention task';
            %'Intensity condition'
            'Expectancy condition'
            %'Group';
            };
        %condition labels, i.e. levels of each condition, in the same order as in
        %the SPM design matrix. One row per factor. second column in plotting order
        P.cval={
            %{'Low','Medium'},[1 2]
            {'Low Exp','Low Exp (post-High)','High Exp'},[1 2]
            %{'Cue 1 Low, Cue 2 Low','Cue 1 High, Cue 2 Low'},[1 2]
            %{'High PCS','Low PCS'},[1 2]
            %{'No Task','Task'},[1 2]
            };

        % clusters to plot
        P.plotclus = {'c1_spm','c2_spm'};
        
     case 2
        %cluster directory name, which also specifies the constrast that will be
        %plotted (i.e. the characters before the underscore)
        P.clusdir='Exp_clusters';
        %factor(s) to plot - if more than one, first factor levels will be on separate plots
        %must be the same characters as used to name the factors in design_batch
        %P.facplot={'Int','Exp'};
        %P.facplot={'Exp','Grp'};
        %P.facplot={'Att','Exp'};
        P.facplot={'Exp'};
        %Full factor names to use for labelling the plots
        P.fact_names = {
            %'Attention task';
            %'Intensity condition'
            'Expectancy condition'
            %'Group';
            };
        %condition labels, i.e. levels of each condition, in the same order as in
        %the SPM design matrix. One row per factor. second column in plotting order
        P.cval={
            %{'Low','Medium'},[1 2]
            {'Low Exp','Low Exp (post-High)','High Exp'},[1 3]
            %{'Cue 1 Low, Cue 2 Low','Cue 1 High, Cue 2 Low'},[1 3]
            %{'High PCS','Low PCS'},[1 2]
            %{'No Task','Task'},[1 2]
            };

        % clusters to plot
        P.plotclus = {'c2_spm'};
        
     case 3
        %cluster directory name, which also specifies the constrast that will be
        %plotted (i.e. the characters before the underscore)
        P.clusdir='Grp_clusters';
        %factor(s) to plot - if more than one, first factor levels will be on separate plots
        %must be the same characters as used to name the factors in design_batch
        %P.facplot={'Int','Exp'};
        %P.facplot={'Exp','Grp'};
        %P.facplot={'Att','Exp'};
        P.facplot={'Grp'};
        %Full factor names to use for labelling the plots
        P.fact_names = {
            %'Attention task';
            %'Intensity condition'
            %'Expectancy condition'
            'Group';
            };
        %condition labels, i.e. levels of each condition, in the same order as in
        %the SPM design matrix. One row per factor. second column in plotting order
        P.cval={
            %{'Low','Medium'},[1 2]
            %{'Low Exp','Low Exp (post-High)','High Exp'},[1 3]
            %{'Cue 1 Low, Cue 2 Low','Cue 1 High, Cue 2 Low'},[1 2]
            {'High PCS','Low PCS'},[1 2]
            %{'No Task','Task'},[1 2]
            };

        % clusters to plot
        P.plotclus = {'c1_spm','c2_spm'};
        
    case 4
        %cluster directory name, which also specifies the constrast that will be
        %plotted (i.e. the characters before the underscore)
        P.clusdir='Grp_Exp_clusters';
        %factor(s) to plot - if more than one, first factor levels will be on separate plots
        %must be the same characters as used to name the factors in design_batch
        %P.facplot={'Int','Exp'};
        P.facplot={'Exp','Grp'};
        %P.facplot={'Att','Exp'};
        %P.facplot={'Grp'};
        %Full factor names to use for labelling the plots
        P.fact_names = {
            %'Attention task';
            %'Intensity condition'
            'Expectancy condition'
            'Group';
            };
        %condition labels, i.e. levels of each condition, in the same order as in
        %the SPM design matrix. One row per factor. second column in plotting order
        P.cval={
            %{'Low','Medium'},[1 2]
            {'Low Exp','Low Exp (post-High)','High Exp'},[1]
            {'High PCS','Low PCS'},[1 2]
            %{'No Task','Task'},[1 2]
            };

        % clusters to plot
        P.plotclus = {'c1_spm','c3_spm'};
        
        % show time windows from all other clusters on this plot?
end
gplot_waveforms(P)
