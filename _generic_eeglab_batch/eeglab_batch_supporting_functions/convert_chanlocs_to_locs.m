function convert_chanlocs_to_locs%(file)

file = 'C:\Data\CORE\eeg\ana\prep\chanlocs.mat'
load(file)

fields  = fieldnames(chanlocs);
keep=zeros(length(fields),1);
keep(1:4)=1;
locs=chanlocs;
locs = rmfield(locs,fields(~keep));
locs = struct2table(locs);

writetable(locs,'chan.txt','Delimiter','\t')