% plot data as dots and boxplots
P.pth = 'C:\Users\cab79\Google Drive\2. Cambridge University\Projects\CRPS ERP study\Winter2017analysis';
% xlsx files must be LONG format. Values in each factor columns should be
% text labels (not numerical) as these will be used to label the charts.
P.filenames = {'Acc_data_longformat.xlsx','RT_data_longformat.xlsx'};
% name of data column
P.datnames = {'Accuracy','ResponseTime'};
% y axis label (should correspond to data type)
P.yaxisnames = {'Accuracy (%)','Response Time (ms)'};
%factor names:
P.xaxisnames = {'Digit','Digit'}; % header of factor represented on x axis
P.subplotnames = {'Side','Side'}; % OPTIONAL factor that splits x axis into subplots
P.groupingnames = {'Group','Group'}; % factor of main interest for comparison - plots as different colours
% other options:
P.groupcolours = [1 0 0; 0 0 1]; % colours of each group on each plot: if two groups, two colours are needed.
P.textsize = 18; % size of text

%filenames = {'N132_data_longformat.xlsx','N124_data_longformat.xlsx','P268_data_longformat.xlsx'};
%datnames = {'N132','N124','P268'};
%yaxisnames = {'N132 Amplitude (µV)','N124 Amplitude (µV)','P268 Amplitude (µV)'};

plotdotbox(P)