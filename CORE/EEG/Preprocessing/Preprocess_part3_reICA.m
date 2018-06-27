

clear all
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
files = dir('*cleaned.set');
cd(filepath);

files_ana = [1]% 1:length(files);

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_loadset('filename',lname,'filepath',filepath);
    
   numcomp = numcompeig(EEG)
   EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',numcomp);
   EEG.reject.gcompreject = ones(1,numcomp);
   sname = [C{1} '_' C{2} '_cleaned_ICA.set'];
   EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 

   clear EEG
end
eeglab
