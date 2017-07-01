breaks = [];
inds = [];
for ev = 1:length(EEG.event)
    br = strcmp(EEG.event(ev).type,'RESP');
    if br
        breaks = [breaks EEG.event(ev).init_time];
        inds = [inds ev];
    end 
end
last = EEG.event(ev).init_time;
etimes = [0 breaks(end)+1 last];