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
D.pdatfile = 'C:\Data\Catastrophising study\Behavioural\All_PR_long_std.xlsx';
% columns to sort INPUT by, in order. negative will sort descending.
D.sortby_in = [1]; 
% factors to sort OUTPUT by, in order (can include old and new factors). Negative will sort descending.
D.sortby_out = [1];
% number of factors in the input (should all be leftmost columns):
D.Nfact_in = 2;
% number of data columns in the input (should all be rightmost column):
D.Ndat_in = 1;
% number of header columns in input data
D.Nhead = 1; 
% OUTPUT factor numbers to average data over
D.avg = []; 
% do std instead of average on the above factors
D.std = 0; 
% average over any non-unique rows in output
D.mean_nonunique=0;

%% enter info for converting from LONG to SHORT (only runs if data is currently long format):
% columns of factors to remove:
D.fact_rm = [1]; 

%% enter info for converting from SHORT to LONG (only runs if data is currently short format):
% name the new factors to create:
D.fact_new = {'S'}; 
% number of levels of each new factor
D.Nlev_fact_new = [53]; 

%% create factor averages for data outputs in short format
% turn on this feature (1 = on, 0 = off)
D.factavg = 0; 
% name the factor averages to create (must be *all* the factors in short form in the output):
D.factavg_new = {}; % strings must be contained within the factor names in the file header
% number of levels of each new factor
D.factavg_lev = []; 

%% apply operation (e.g. subtraction) to levels of one factor for data outputs in short format
% 1 = turn on, 0 = turn off
D.applyop = 0;
% name ALL factors in short format in the output:
D.factshort = {'Change Probability','Side','Digit Change'}; % strings must be contained within the factor names in the file header
% number of levels of each current short-form factor
D.factshort_lev = [3 2 2]; 
% index of above factor(s) to apply operation to
D.factop = 3;
% define operation
D.op = 'i1 - i2'; % subtraction of level 1 minus level 2
%D.op = '(i2 - i1) ./ i2'; % percentage change
% operation name
D.opname = 'subtract_DC';


%% run the function
convertLS(D)