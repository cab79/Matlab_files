addpath(genpath('M:\Matlab\Matlab_files\SA study'));
filepath = 'C:\Data\Sustained attention study\Results'; 
subpath = '\EEG data\Part 2';
savepath = 'C:\Data\Sustained attention study\Preprocessed';
parts = [1];
combine=0; 

subj = dir(fullfile(filepath,'S*'));
for s = 1:length(subj)
    dataimport_cnt(fullfile(filepath,subj(s).name,subpath),'2.cnt',parts,combine,savepath,subj(s).name)
end