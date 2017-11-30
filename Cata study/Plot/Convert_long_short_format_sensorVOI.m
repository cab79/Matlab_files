% CONVRSION FROM LONG TO SHORT FORMAT AND VICE VERSA
% "Long" format data has columns for different factors in the leftmost columns,
% and a single rightmost column containing the data.
% "Short" format data can have some or no factor columns, and has a number of rightmost columns
% containing data that are combinations of factors. 

% To convert between 2 different "Short" formats, the script needs to be run twice:
% 1. Short to Long
% 2. Long to new Short

clear all
%% file info
D.format = 'long'; %CURRENT format of data
% file path and name
D.pdatfile = 'C:\Data\Catastrophising study\SPMstats\Include1\Between-SPN\t-3000_0_b-3000_-2500_Grp_Att_Subject_orig_cleaned_SPN_spm\Grp_clusters\VOI_c9.xlsx';
% columns to sort INPUT by, in order. negative will sort descending.
D.sortby_in = []; 
% factors to sort OUTPUT by, in order (can include old and new factors). Negative will sort descending.
D.sortby_out = [];
% number of factors in the input (should all be leftmost columns):
D.Nfact_in = 3;
% number of data columns in the input (should all be rightmost column):
D.Ndat_in = 1;
% number of header rows in input data
D.Nhead = 1; 
% OUTPUT factor numbers to average data over
D.avg = [3]; 
% average over any non-unique rows in output
D.mean_nonunique=0;

%% enter info for converting from LONG to SHORT (only runs if data is currently long format):
% columns of factors to remove:
D.fact_rm = []; 

%% enter info for converting from SHORT to LONG (only runs if data is currently short format):
% name the new factors to create:
D.fact_new = {'Task','Stim','Exp'}; 
% number of levels of each new factor
D.Nlev_fact_new = [2 2 2]; 

%% run the function
convertLS(D)