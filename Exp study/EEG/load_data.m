addpath(genpath('M:\Matlab\Matlab_files\Exp study'));
filepath = 'C:\Data\Expectancy Study\raw data'; 
subpath = '';
savepath = 'C:\Data\Expectancy Study\Preprocessed';
%parts = [1];
combine=1; 

subj = dir(fullfile(filepath,'*'));
for s = 11%3:length(subj)
    parts = length(dir(fullfile(filepath,subj(s).name,subpath,'Block*.cnt')));
    dataimport_cnt(fullfile(filepath,subj(s).name,subpath),'Block*.cnt',parts,combine,savepath,subj(s).name)
end