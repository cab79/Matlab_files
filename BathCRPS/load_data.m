addpath(genpath('M:\Matlab\Matlab_files\BathCRPS'));
filepath = 'C:\Users\Chris\Google Drive\1. Liverpool University\Projects\BathCRPS'; 
savepath = 'C:\Users\Chris\Google Drive\1. Liverpool University\Projects\BathCRPS';
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'*00171a_ERP Artifact Rejection.vhdr'));
for s = 1:length(subj)
    dataimport_bv(filepath,subj(s).name,parts,combine,savepath)
end