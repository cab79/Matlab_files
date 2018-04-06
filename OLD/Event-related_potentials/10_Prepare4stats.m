clear all

chan_locs = 'C:\Documents and Settings\mdmoscab\Desktop\Data analysis files\EEG analysis Matlab\chan_MR62.locs';
[tmpeloc labels Th Rd indices] = readlocs(chan_locs,'filetype','loc');
labels = upper(labels);
p1 = find(ismember(labels,{'AF7','F3','F5','F7'}));
p2 = find(ismember(labels,{'FZ','F1','F3','AFZ','AF1','AF3','FPZ','FP1','FP3'}));
p3 = find(ismember(labels,{'AF8','F4','F6','F8'}));
p4 = find(ismember(labels,{'FT7','FC5','FC3','T7','C5','C3','TP7','CP5','CP3'}));
p5 = find(ismember(labels,{'FC1','FC2','CZ','C1','C2','CPZ','CP1','CP2'}));
p6 = find(ismember(labels,{'FT8','FC6','FC4','T8','C6','C4','TP8','CP6','CP4'}));
p7 = find(ismember(labels,{'P7','P3','P5','PO7'}));
p8 = find(ismember(labels,{'PZ','P1','P2','POZ','PO3','PO4','OZ','O1','O2'}));
p9 = find(ismember(labels,{'P8','P4','P6','PO8'}));

files = {'gamma1.mat';
'gamma2.mat';
'beta.mat'; 
'alpha.mat'; 
'theta.mat'; 
'delta.mat'}; 

for f = 1:length(files)
    
load(files{f});
[pth,fname,ext] = fileparts(files{f});
eval(['data = ' fname ';'])

for d = 1:9
    eval(['ndata(' num2str(d) ',:,:,:) = mean(data(p' num2str(d) ',:,:,:),1)']);
end

eval(['save ' fname '_9reg.mat ndata']);

end
