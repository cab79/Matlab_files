clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*ICA.set');


files_ana = 1:length(files);
comtab=NaN(length(files_ana),60);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    comtab(f,1:length(EEG.reject.gcompreject))=EEG.reject.gcompreject;
end

nrej = nansum(comtab,2);
outfiles = files(find(nrej>=10))