addpath(genpath('M:\Matlab\Matlab_files\BathCRPS'));
filepath = 'C:\Data\BathCRPS'; 
savepath = 'C:\Data\BathCRPS';
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'*Average.vhdr'));
for s = 1:length(subj)
    dataimport_bv(filepath,subj(s).name,parts,combine,savepath)
end