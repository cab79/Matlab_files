%addpath(genpath('Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory')); 
dbstop if error
filepath = 'C:\Data\NTIP\Raw'; %where the original data is stored
savepath = 'C:\Data\NTIP\Raw'; %where the sets will go
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'*.vhdr'));
for s = 1:length(subj)
    dataimport_bv(filepath,subj(s).name,parts,combine,savepath)
end