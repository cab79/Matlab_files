function Fig4_SensorWaveforms
%% plot ERP waveforms from clusters using the amazing Gramm toolbox!

clear all
close all

% save options
P.save_waveforms=1;
P.save_topo=1;

%% figures to plot
P.plotnum = 4;

for pnum = P.plotnum

    switch pnum
        case 1
            P=ERP_settings;
            P(p).xplotpos = 1;
            P(p).yplotpos = 1;
            
            %cluster directory name, which also specifies the constrast that will be
            %plotted (i.e. the characters before the underscore)
            P.clusdir='ExpAC_clusters';
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
                'Cue'
                %'Group';
                };
            %condition labels, i.e. levels of each condition, in the same order as in
            %the SPM design matrix. One row per factor. second column in plotting order
            P.cval={
                %{'Low','Medium'},[1 2]
                {'Low','Low (prior High)','High'},[1 3]
                %{'Cue 1 Low, Cue 2 Low','Cue 1 High, Cue 2 Low'},[1 3]
                %{'High PCS','Low PCS'},[1 2]
                %{'No Task','Task'},[1 2]
                };

            % clusters to plot
            P.plotclus = {'c1_spm','c2_spm','c3_spm','c4_spm'};
            
            gplot_waveforms(P)
            
            
            %-------- plot time series -----------%
            P.gtitle = clnames{cl};
            g = gplot_timeseries(P)
            %-------------------------------------%


         case 2
            P=ERP_settings;
             
            %cluster directory name, which also specifies the constrast that will be
            %plotted (i.e. the characters before the underscore)
            P.clusdir='ExpACGrp_clusters';
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
                {'High PC','Low PC'},[1 2]
                %{'No Task','Task'},[1 2]
                };

            % clusters to plot
            P.plotclus = {'c1_spm','c2_spm'};
            
            gplot_waveforms(P)


        case 3
            P=ERP_settings;
            
            %cluster directory name, which also specifies the constrast that will be
            %plotted (i.e. the characters before the underscore)
            P.clusdir='Grp_ExpAC_clusters';
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
                'Cue'
                'Group';
                };
            %condition labels, i.e. levels of each condition, in the same order as in
            %the SPM design matrix. One row per factor. second column in plotting order
            P.cval={
                %{'Low','Medium'},[1 2]
                {'Low','Low (prior High)','High'},[1]
                {'High PC','Low PC'},[1 2]
                %{'No Task','Task'},[1 2]
                };

            % clusters to plot
            P.plotclus = {'c1_spm'};

            gplot_waveforms(P)

        case 4
            P=ERP_settings;
            
            %cluster directory name, which also specifies the constrast that will be
            %plotted (i.e. the characters before the underscore)
            P.clusdir='Grp_ExpAC_clusters';
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
                'Cue'
                'Group';
                };
            %condition labels, i.e. levels of each condition, in the same order as in
            %the SPM design matrix. One row per factor. second column in plotting order
            P.cval={
                %{'Low','Medium'},[1 2]
                {'Low','Low (prior High)','High'},[3]
                {'High PC','Low PC'},[1 2]
                %{'No Task','Task'},[1 2]
                };

            % clusters to plot
            P.plotclus = {'c1_spm'};

            gplot_waveforms(P)

    end
end
end

function P=ERP_settings
%% generic directories for all analyses for this study
%-------------------------------------------------------------
% directory in which SPM is saved
%P.spm_path = 'C:\Data\Catastrophising study\SPMstats\Source\1_grp\_Time_Int_Exp_Subject_spm_t416_478';
%P.wavetype = 'source'; % source or sensor?
P.spm_path = 'C:\Data\Catastrophising study\SPMstats\Include1\Between\t-3000_-2_b-3000_-2500_m_-2500_-1000_Grp_Exp_Subject_orig_cleaned_SPNall_spm';
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
P.xaxisname = {'post-cue time (ms)'};
P.yaxisname = {'amplitude (uV)'};
% time in ms
P.xlimits = [-2700 -1000];
% vertical dashed lines to indicate events, time in ms
P.xlinedashed = [-5000 -2500 0];
% polygon times (or leave blank to extract from clustertable)
P.poly = [];
% change zero latency to this time
P.timezero = -2500;

% plot topographies
P.plot_topo=1;
P.eventtypes = {'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'};
P.eeglab_path = 'C:\Data\Catastrophising study\Preprocessed';
P.no_plot_ele=[];
P.topo_subtimewin=0;%[-2265 -2265]; % time window length. single value: plots multiple topos in consecutive windows with this length. 2 values: specifies a window. 0 = whole cluster mean.
end