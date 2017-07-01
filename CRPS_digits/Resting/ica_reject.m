clear all
loadpaths
typerange = {'RHAN','LHAN','BELB','RELX'};
files = dir(fullfile(filepath,'*100Hz.Exp3.set'));
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab 

for f = 1:length(files)

    filename = files(f).name;
    [pth nme ext] = fileparts(filename);
    for e = 1:length(typerange)
        EEG = pop_loadset('filename',[nme '_' typerange{e} ext],'filepath',filepath);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); % copy it to ALLEEG
        EEG = pop_rejepoch(EEG, [], 0);
        if e==1; EEG = pop_selectcomps(EEG, [1:35]); pause; end
        rcomp = EEG.reject.gcompreject;
        if ~exist('rejcomp','var'); rejcomp = rcomp; end
        rejcomp = (rcomp | rejcomp);
        clear EEG
    end
    
    for e = 1:length(typerange)
        EEG = ALLEEG(e);
        CURRENTSET = e;
        EEG.reject.gcompreject = rejcomp;
        EEG = pop_subcomp(EEG,find(rejcomp),0);
    end
    eeglab redraw
    keyboard;
    
     for e = 1:length(typerange)
        EEG = ALLEEG(e);
        CURRENTSET = e;
        EEG = pop_rejepoch(EEG);
        %pop_eegplot(EEG,1,1,1);
        
        pop_saveset(EEG,'filename',[nme '_' typerange{e} '_subcomp' ext],'filepath',filepath);
        %clear EEG;
    end
    
    clear EEG rejcomp rcomp
    [ALLEEG] = pop_delset(ALLEEG,1:length(ALLEEG)); 
end

