%% plot VOI data as dots and boxplots

% PREREQUISITES: 
% - VOI data file(s) in excel LONG format (e.g. output from Convert_VOImat_to_excel.m)
close all
P.pth = 'C:\Data\CORE\SPMstats\t-200_299_b-200_0_m_0_299_CP_Odd_DC_Subject_4_cleaned_tm_spm\Odd_clusters';
%P.pth = 'C:\Data\CORE\SPMstats\t-200_299_b-200_0_m_0_299_CP_Grp_Odd_Subject_4_cleaned_tm_spm\Odd_clusters';

% xlsx files must be LONG format. Values in each factor columns should be
% text labels (not numerical) as these will be used to label the charts.
P.filenames = {'VOI_c1_long_20170804T221045.xlsx'}; % use * to load multiple files, each with separate plots. Filename must have extension
%P.filenames = {'VOI_c1_long_20170804T192353.xlsx'}; % use * to load multiple files, each with separate plots. Filename must have extension
% name of data column
P.datnames = {'Data'};
% y axis label (should correspond to data type)
P.yaxisnames = {'mean cluster activity (µV)'};
%factor names:
%P.xaxisnames = {'Stimulus Intensity'}; % header of factor represented on x axis
P.xaxisnames = {'Digit Change'}; % header of factor represented on x axis
%P.subplotnames = {}; % OPTIONAL factor that splits x axis into subplots
P.subplotnames = {'Change Probability'}; % OPTIONAL factor that splits x axis into subplots
    P.sublevel = 0.1;% OPTIONAL selection of level with subplotnames factor - necessary for line plots
%P.groupingnames = {'Expectation Cues'}; % factor of main interest for comparison - plots as different colours
P.groupingnames = {'Oddball effect'}; % factor of main interest for comparison - plots as different colours
%P.groupingnames = {'Group'}; % factor of main interest for comparison - plots as different colours
% other options:
%colours = [0 0.8 1; 1 0.2 0]; % blue, red
P.groupcolours = [0.8 0.2 0.8; 0.2 1 0.2]; % green, purple 
%P.groupcolours = []; % colours of each group on each plot: if two groups, two colours are needed.
P.textsize = 20; % size of text
P.label_points = 0; % 1 = labels points with subject name; 0 = no labelling
P.plot_sep_fg=0;

% 1 = plot, 0 = no plot. Will align plots horizonally in figure, so a small number should be slected (2-3).
P.plottypes = {
    'jitter', 0;
    'boxplot', 1;
    'violin', 0;
    'confidence',0;
    'line',0;
    };
    
P.save_figure=1;

%% call function
plotdotbox(P)