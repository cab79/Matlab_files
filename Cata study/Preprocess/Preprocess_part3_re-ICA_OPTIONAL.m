clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*cleaned.set');
addpath(genpath('M:\Matlab\Matlab_files\Cata study'));
%files_ana = [38,39];%usable if further cleaned
%files_ana = [12,17,23,28,35];%warnings or unusable
%files_ana = [1:11 13:16 18:20 24 26 29 31:34 36:37 40]; % the rest
files_ana = [14:16 18:20 24 26 29 31:34 36:37 40]; % the rest temp

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_loadset('filename',lname,'filepath',filepath);
    
   EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',30);
   sname = [C{1} '_' C{2} '_cleaned_ICA.set'];
   EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

   clear EEG
end
eeglab
