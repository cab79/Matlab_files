clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*23_orig_cleaned.set');
load('C:\Data\Catastrophising study\Orig\chanlocs.mat');

basebin = [-5.5 -5];
stimtypes = {'c0','c1','c2','c3','c4','c5','c6','c7','c8'};

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    EEG = FTrejman(EEG,[0 0]);
    
end