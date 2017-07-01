clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*ICA.set');

comtab=[];
files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    comtab(f,:)=EEG.reject.gcompreject;
end

nrej = sum(comtab,2);
outfiles = files(find(nrej>=10))