clear all
delay = 50; % ms - e.g. adds this to all conditions. caused by delay in EEG marker relative to stimuli
jitter = 180; %ms - removes this from the start of the EEG from the below conditions
cut_start = -3000; %ms
paste_end = -2;
etypes = {'c2b','c4b','c6b','c8b'};

filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*cleaned_SPNpn.set');

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    % fix time delay
    dpdelay = delay*(EEG.srate/1000);
    if delay>0
        EEG.data = cat(2,EEG.data(:,1:dpdelay,:), EEG.data(:,1:end-dpdelay,:));
    else
        EEG.data = cat(2,EEG.data(:,-dpdelay+1:end,:), EEG.data(end-dpdelay+1:end,:));
    end
    
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