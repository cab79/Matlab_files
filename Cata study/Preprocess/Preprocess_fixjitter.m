clear all
delay = 20; % ms - e.g. caused by delay in EEG marker relative to stimuli
jitter = 242; %ms
cut_start = -5500; %ms
paste_end = -2;
etypes = {'c1','c3','c5','c7'};

filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*cleaned*.set');

files_ana = [39 40] %2:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    % fix time delay
    dpdelay = delay*(EEG.srate/1000);
    EEG.data = cat(2,EEG.data(:,1:dpdelay,:), EEG.data(:,1:end-dpdelay,:));
    
    % fix jitter
    dpjitter = jitter*(EEG.srate/1000);
    dp_cut_start = find(EEG.times==cut_start);
    dp_paste_end = find(EEG.times==paste_end);
    
    for ev = 1:length(EEG.event)
        if ismember(EEG.event(ev).type,etypes)
            cut = EEG.data(:,dp_cut_start:dp_cut_start+dpjitter-1,ev);
            bc_val = mean(EEG.data(:,dp_paste_end-dpjitter,ev),2) - mean(cut,2);
            bc_val = repmat(bc_val,1,size(cut,2));
            EEG.data(:,:,ev) = cat(2,EEG.data(:,dp_cut_start+dpjitter:dp_paste_end,ev), cut+bc_val, EEG.data(:,dp_paste_end+1:end,ev));
        end
    end
    
    EEG = pop_saveset(EEG,'filename',files(f).name,'filepath',filepath); 
end