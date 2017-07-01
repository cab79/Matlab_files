addpath(genpath('m:\Matlab\Matlab_files\Cata study'));
filepath = 'C:\Data\Catastrophising study\BrainvisionData\Study2EEGraw'; 
%subpath = '\EEG data\Part 2';
savepath = 'C:\Data\Catastrophising study\Orig\';
parts = [1];
combine=1; 

group = [1]
for s = 1:length(group)
    subj = num2str(group(s));
    if length(subj)<2
        subj = ['0' subj]
    end
    subjname = ['C' subj '*'];
    savename = ['C' subj];
    dataimport_cnt(filepath,subjname,parts,combine,savepath,savename)
end