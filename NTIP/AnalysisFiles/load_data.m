addpath(genpath('Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory')); 
filepath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Original Data'; %where the original data is stored
savepath = 'Y:\Marie Shorrock\NTIP\Pilot_Tim_Auditory\Set'; %where the sets will go
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'NTIP_0000_Tim*.vhdr'));
for s = 1:length(subj)
    dataimport_bv(filepath,subj(s).name,parts,combine,savepath)
end