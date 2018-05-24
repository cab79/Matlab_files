addpath(genpath('C:\Data\Matlab\Matlab_files\SCAP\draft_import'));
filepath = 'C:\Data\SCAP\EEG\Raw'; 
savepath = 'C:\Data\SCAP\EEG\Raw';
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'0021a_DC Detrend.vhdr'));
for s = 1:length(subj)
    dataimport_bv(filepath,subj(s).name,parts,combine,savepath)
end