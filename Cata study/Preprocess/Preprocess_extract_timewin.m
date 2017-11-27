% Extracts a timewindow from epoched data that si specific to trial type
% Saves output as a new epoched data file.

clear all

curr_file = '*cleaned.set';
new_file_ext = '_SPN';
% etypes:
% col 1: trial type
% col 2: timewin in current data
% col 3: new definition of trial start/end
etypes = {
    {'c0','c1','c3','c5','c7'},{[-3000, 0]},{[-3000, 0]};
    {'c2','c4','c6','c8'},{[-5500, -2500]},{[-3000, 0]}
    };

filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir(curr_file);

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    olddata = EEG.data;
    newdata=[];
    
    for et = 1:size(etypes,1)
    
        old_tw = etypes{et,2}{:}; %ms
        new_tw = etypes{et,3}{:}; %ms

        dp_cut_start = find(EEG.times==old_tw(1));
        dp_cut_end = find(EEG.times==old_tw(2));
        new_start = find(EEG.times==new_tw(1));
        new_end = find(EEG.times==new_tw(2));

        for ev = 1:length(EEG.event)
            if ismember(EEG.event(ev).type,etypes{et,1})
                newdata(:,:,ev) = olddata(:,dp_cut_start:dp_cut_end-1,ev);
            end
        end
    end
    
    EEG.data = newdata;
    EEG.pnts = size(EEG.data,2);
    EEG.xmin = new_tw(1)/1000;
    EEG.xmax = new_tw(2)/1000 - 1/EEG.srate;
    EEG.times = EEG.times(new_start:new_end-1);
    
    [~,savename,~] = fileparts(files(f).name);
    savename = [savename new_file_ext];
    EEG = pop_saveset(EEG,'filename',savename,'filepath',filepath); 
end