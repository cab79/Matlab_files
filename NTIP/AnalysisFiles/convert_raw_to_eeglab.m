% This simple script identifies files for conversion to EEGLAB format for a
% specific study and sets some parameters for the conversion. It calls the generic function
% "import_eeglab" to do the actual conversion.

%% SETUP
dbstop if error % optional instruction to stop at a breakpoint if there is an error - useful for debugging
S.loadpath = 'C:\Users\cab79\Desktop'; %where the raw input data is stored
S.savepath = 'C:\Users\cab79\Desktop'; %where the .set output files will go
S.inputfileext = 'cnt'; % file extension of input data: supports 'vhdr' (Brainvision), 'cnt' (Neuroscan)
S.outputfileprefix = ''; % suffix to add to output file, if needed
S.outputfilesuffix = ''; % suffix to add to output file, if needed
S.chanexcl = [31,32]; % exclude channels, or leave empty as []
S.addchanloc = 'C:\Data\Matlab\eeglab13_6_5b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp'; % add channel locations from this path; or leave as ''

%% RUN
import_eeglab(S)