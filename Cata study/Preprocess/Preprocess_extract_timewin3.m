% Extracts a timewindow from epoched data that si specific to trial type
% Saves output as a new epoched data file.

clear all

curr_file = '*cleaned.set';
new_file_ext = '_SPNall';
% etypes:
% col 1: trial type
% col 2: timewin in current data
% col 3: new definition of trial start/end
% col 4: append to event name
etypes = {
    {'c0','c1','c3','c5','c7'},{[-3000, 0]},{[-3000, 0]},'';
    {'c2','c4','c6','c8'},{[-3000, 0]},{[-3000, 0]},'a';
    {'c2','c4','c6','c8'},{[-5500, -2500]},{[-3000, 0]},'b'
    };

filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir(curr_file);

files_ana = 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEGall = pop_loadset('filename',files(f).name,'filepath',filepath);
    
    for et = 1:size(etypes,1)
    
        EEG(et) = pop_selectevent(EEGall,'type',etypes{et,1});

        olddata = EEG(et).data;
        newdata=[];
    
        old_tw = etypes{et,2}{:}; %ms
        new_tw = etypes{et,3}{:}; %ms

        dp_cut_start = find(EEG(et).times==old_tw(1));
        dp_cut_end = find(EEG(et).times==old_tw(2));

        for ev = 1:length(EEG(et).event)
            if ismember(EEG(et).event(ev).type,etypes{et,1})
                newdata(:,:,ev) = olddata(:,dp_cut_start:dp_cut_end-1,ev);
                EEG(et).epoch(ev).eventtype = [EEG(et).epoch(ev).eventtype etypes{et,4}];
                EEG(et).event(ev).type = [EEG(et).event(ev).type etypes{et,4}];
                %EEG(et).event(ev).latency = EEG(et).event(ev).latency+(et-1);
            end
        end
        
        EEG(et).data = newdata;
        EEG(et).pnts = size(EEG(et).data,2);
        EEG(et).xmin = new_tw(1)/1000;
        EEG(et).xmax = new_tw(2)/1000 - 1/EEG(et).srate;
        EEG(et).times = EEG(et).times(dp_cut_start:dp_cut_end-1);
        
        % merge data structures
        if et==1
            EEGout = EEG(et);
        else
            EEGout.data = cat(3,EEGout.data,EEG(et).data);
            EEGout.trials = EEGout.trials+EEG(et).trials;
            EEGout.epoch = [EEGout.epoch,EEG(et).epoch];
            EEGout.event = [EEGout.event,EEG(et).event];
            EEGout.urevent = [EEGout.urevent,EEG(et).urevent];
            
            t=num2cell(1:length([EEGout.epoch.event])); [EEGout.epoch.event]=t{:};
            t=num2cell(1:length([EEGout.event.epoch])); [EEGout.event.epoch]=t{:};
            
        end
    end
    
    
    [~,savename,~] = fileparts(files(f).name);
    savename = [savename new_file_ext];
    pop_saveset(EEGout,'filename',savename,'filepath',filepath); 
end