clear all
delay = 0; % ms - e.g. adds this to all conditions. caused by delay in EEG marker relative to stimuli
jitter = [166, -84, 100]; %ms - removes this from the start of the EEG from the below conditions
cut_start = -3000; %ms
paste_end = -2;
etypes = {'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'};
use_etype = [1 2 3 4; 5 6 7 8; 9 10 11 12];

filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*cleaned_SPNall.set');

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    % fix time delay
    dpdelay = delay*(EEG.srate/1000);
    if delay>0
        EEG.data = cat(2,EEG.data(:,1:dpdelay,:), EEG.data(:,1:end-dpdelay,:)); % shift data to right
    elseif delay<0
        EEG.data = cat(2,EEG.data(:,-dpdelay+1:end,:), EEG.data(:,end+dpdelay+1:end,:)); % shift data to left
    end
    
    % fix jitter
    dpjitter = jitter*(EEG.srate/1000);
    dp_cut_start = find(EEG.times==cut_start);
    dp_paste_end = find(EEG.times==paste_end);
    
    for ev = 1:length(EEG.event)
        for et = 1:size(use_etype,1)
            if ismember(EEG.event(ev).type,etypes(use_etype(et,:)))
                if dpjitter(et)>0 % shift left
                    cut = EEG.data(:,dp_cut_start:dp_cut_start+dpjitter(et)-1,ev);
                    bc_val = mean(EEG.data(:,dp_paste_end-dpjitter(et),ev),2) - mean(cut,2);
                    bc_val = repmat(bc_val,1,size(cut,2));
                    EEG.data(:,:,ev) = cat(2,EEG.data(:,dp_cut_start+dpjitter(et):dp_paste_end,ev), cut+bc_val, EEG.data(:,dp_paste_end+1:end,ev));
                else % shift right
                    cut = EEG.data(:,dp_paste_end+dpjitter(et)+1:dp_paste_end,ev);
                    bc_val = mean(EEG.data(:,dp_cut_start,ev),2) - mean(cut,2);
                    bc_val = repmat(bc_val,1,size(cut,2));
                    EEG.data(:,:,ev) = cat(2,EEG.data(:,1:dp_cut_start,ev),cut+bc_val,EEG.data(:,dp_cut_start+1:dp_paste_end+dpjitter(et),ev), EEG.data(:,dp_paste_end+1:end,ev));
                end
            end
        end
    end
    
    EEG = pop_saveset(EEG,'filename',files(f).name,'filepath',filepath); 
    %EEG = pop_saveset(EEG,'filename',['temp_' files(f).name],'filepath',filepath); 
end