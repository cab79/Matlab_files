addpath(genpath('M:\Matlab\Matlab_files\PET-LEP'));
filepath = 'C:\Data\PET-LEP\RawData'; 
savepath = 'C:\Data\PET-LEP\Preprocessed';
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'*part2.vhdr'));
for s = 1:length(subj)
    dataimport_bv(filepath,subj(s).name,parts,combine,savepath)
end